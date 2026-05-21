%% Retinal microglia: Bayesian posterior predictive bands (eta_O fixed)
% Model: microglia_model_opt12 -> [OPL, IPL, GCL, R, B]
% Free params (positive): [eta_G, k_OI, k_IO, k_IG, k_G], plus noise sigma
% Prior: flat on log-parameters within very wide bounds (no biology priors)
% Noise model: y = model + N(0, sigma^2)
% Output: 5–95% posterior predictive bands + 95% credible intervals
clear; clc; rng(1);

%% ---------------------- Load data ----------------------
load('Retinal_data.mat');   % time, OPL_individual, IPL_individual, GCL_NFL_individual
tdata    = time(:);
OPL_data = OPL_individual;
IPL_data = IPL_individual;
GCL_data = GCL_NFL_individual;

rowmean  = @(A) mean(A,2,'omitnan');
OPL_mean = rowmean(OPL_data);
IPL_mean = rowmean(IPL_data);
GCL_mean = rowmean(GCL_data);

%% ---------------------- Fixed params & ICs ----------------------
p = struct();
p.k_R   = 1.71;
p.F     = 0.03;
p.Vol_P = 0.06;
p.Vol_B = 0.9548;
p.k_B   = 0.0578;
p.eta_O = 1667;        % fixed

ICs   = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2021000, 0];
tspan = [tdata(1), tdata(end)];

% NEW: initial G value for logistic term
p.G0  = ICs(3);   % or whatever you want to use as G(0)

%% ---------------------- Build observation vectors ----------------------
[obsVec, timeVec, layerVec] = pack_observations(OPL_data, IPL_data, GCL_data, tdata);
Nobs = numel(obsVec);

%% ---------------------- LS fit (for good starting point) ----------------------
q0  = [266, 0.264811, 0.488318, 0.112699, 0.0734471];  % [eta_G, k_OI, k_IO, k_IG, k_G]
lb  = [1e-8, 1e-8, 1e-8, 1e-8, 1e-8];
ub  = [1e8,  1e2,  1e2,  1e2,  1e2];

opts = optimoptions('lsqnonlin','Display','iter','FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',2e4,'MaxIterations',800);

warnID = 'MATLAB:ode15s:IntegrationTolNotMet';
st = warning('query', warnID); warning('error', warnID);

obj = @(q) residuals_vec(q, obsVec, timeVec, layerVec, ICs, tspan, p);
[q_hat, ~, ~, ~, ~, ~, J] = lsqnonlin(obj, q0, lb, ub, opts);
r_hat  = obj(q_hat);
SSE    = sum(r_hat.^2);
k_free = numel(q_hat);
sigma2_hat = SSE / max(Nobs - k_free, 1);
sigma_hat  = sqrt(sigma2_hat);
warning(st.state, warnID);

fprintf('\n=== LS (for initialization) with eta_O fixed ===\n');
names = {'eta_G','k_OI','k_IO','k_IG','k_G'};
for i=1:5, fprintf('%-5s = %.6g\n', names{i}, q_hat(i)); end
fprintf('sigma  = %.4g\n', sigma_hat);

%% ---------------------- MCMC settings ----------------------
% We sample theta = log(params) and s = log(sigma)
theta0 = log(q_hat(:));
s0     = log(sigma_hat);

% Very wide "no prior info" bounds (on original scale) -> used as support
BETA_LO = [1e-12, 1e-8, 1e-8, 1e-8, 1e-8];   % lower bounds
BETA_HI = [1e10,  1e3,  1e3,  1e3,  1e3];    % upper bounds

% Proposal covariance for log-parameters: use LS Hessian approx if available
% C_q ≈ sigma^2 * inv(J'J); convert to log-scale via diag(1/q_hat)
Cq = eye(5);
if ~isempty(J)
    H = J.'*J;

    % Symmetrize and ensure FULL for eig (sparse -> full to allow 2 outputs)
    symH = (H + H.')/2;
    if issparse(symH), symH = full(symH); end

    % Eigen-regularize
    [V,D] = eig(symH);
    d = diag(D);
    d = max(d, 1e-10);                         % floor tiny/neg eigenvalues
    Hreg = V*diag(d)*V.';

    % Linear-scale covariance and log-transform
    C_lin = sigma2_hat * (Hreg \ eye(size(Hreg,1)));
    S = diag(1./max(q_hat(:),1e-12));          % d log(q)/dq ≈ 1/q_hat
    Cq = S * C_lin * S;
end

% Scale proposals a bit (tune if needed)
scale = 2.4^2 / numel(theta0);
PropCov = scale * (Cq + 1e-4*eye(5));  % add jitter
PropChol = chol(PropCov, 'lower');
propStd_s = 0.15;                       % step size for log(sigma)

% Chain lengths
n_iter  = 16000;
burn_in = 4000;
thin    = 4;

%% ---------------------- Log-posterior (up to additive const) ----------------------
logpost = @(theta,s) logpost_theta(theta, s, obsVec, timeVec, layerVec, ICs, tspan, p, ...
                                   BETA_LO, BETA_HI);

% Pre-allocate
keep_idx = 0;
n_keep   = floor((n_iter - burn_in)/thin);
TH_SAMPLES = zeros(n_keep, 5);
SIG_SAMPLES= zeros(n_keep, 1);

% Init
theta = theta0; s = s0;
lp_curr = logpost(theta, s);

accept = 0; acc_s = 0;

fprintf('\n=== MCMC (RW-MH) starting ===\n');
for it = 1:n_iter
    % propose log-parameters
    theta_prop = theta + PropChol*randn(5,1);
    s_prop     = s + propStd_s*randn;
    lp_prop    = logpost(theta_prop, s_prop);

    if isfinite(lp_prop)
        a = lp_prop - lp_curr;
        if log(rand) < a
            theta = theta_prop; s = s_prop; lp_curr = lp_prop;
            accept = accept + 1;
        end
    end

    % save after burn-in & thinning
    if it > burn_in && mod(it - burn_in, thin)==0
        keep_idx = keep_idx + 1;
        TH_SAMPLES(keep_idx,:) = theta(:).';
        SIG_SAMPLES(keep_idx)  = s;
    end

    % (optional progress)
    if mod(it, 2000)==0
        fprintf(' iter %5d / %5d, acc rate ~ %.2f\n', it, n_iter, accept/it);
    end
end
acc_rate = accept/n_iter;
fprintf('MCMC done. Acceptance ~ %.2f\n', acc_rate);

% Convert to original scale
Qpost   = exp(TH_SAMPLES);        % [n_keep x 5]
Sigpost = exp(SIG_SAMPLES);       % [n_keep x 1]

%% ---------------------- 95%% credible intervals (parameters) ----------------------
prc = @(x,p) prctile(x,p,1);  % along rows
q_med = prc(Qpost,50);
q_lo  = prc(Qpost,2.5);
q_hi  = prc(Qpost,97.5);
sigma_med = median(Sigpost);
sigma_lo  = prctile(Sigpost,2.5);
sigma_hi  = prctile(Sigpost,97.5);

fprintf('\n=== 95%% credible intervals (posterior) ===\n');
for i=1:5
    fprintf('  %-5s : median=%.6g   95%% CI = [%.6g, %.6g]\n', ...
        names{i}, q_med(i), q_lo(i), q_hi(i));
end
fprintf('  sigma : median=%.6g   95%% CI = [%.6g, %.6g]\n', ...
    sigma_med, sigma_lo, sigma_hi);

%% ---------------------- Posterior predictive bands ----------------------
nsim = min(250000, size(Qpost,1));           % was 600
idx  = randi(size(Qpost,1), nsim, 1);      % sample with replacement (more iid-like)
Qsim = Qpost(idx,:);

tgrid   = linspace(tspan(1), tspan(2), 1000);
opts_ode= odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',(tspan(2)-tspan(1))/500);

Yall_O = zeros(numel(tgrid), nsim);
Yall_I = zeros(numel(tgrid), nsim);
Yall_G = zeros(numel(tgrid), nsim);
Rall   = zeros(numel(tgrid), nsim);
Ball   = zeros(numel(tgrid), nsim);

for k = 1:nsim
    p_k = overwrite_from_q(p, Qsim(k,:));
    [~, Yk] = ode15s(@(tt,yy) microglia_model_opt12(tt,yy,p_k), tgrid, ICs, opts_ode);
    Yall_O(:,k) = Yk(:,1);
    Yall_I(:,k) = Yk(:,2);
    Yall_G(:,k) = Yk(:,3);
    Rall(:,k)   = Yk(:,4);
    Ball(:,k)   = Yk(:,5);
end

% 2.5/50/97.5% bands
qO_lo = prctile(Yall_O,  2.5, 2);  qO_md = prctile(Yall_O, 50, 2);  qO_hi = prctile(Yall_O, 97.5, 2);
qI_lo = prctile(Yall_I,  2.5, 2);  qI_md = prctile(Yall_I, 50, 2);  qI_hi = prctile(Yall_I, 97.5, 2);
qG_lo = prctile(Yall_G,  2.5, 2);  qG_md = prctile(Yall_G, 50, 2);  qG_hi = prctile(Yall_G, 97.5, 2);
qR_lo = prctile(Rall,    2.5, 2);  qR_md = prctile(Rall,   50, 2);  qR_hi = prctile(Rall,   97.5, 2);
qB_lo = prctile(Ball,    2.5, 2);  qB_md = prctile(Ball,   50, 2);  qB_hi = prctile(Ball,   97.5, 2);

%% ---------------------- Plot (poster layout) ----------------------
%% ---------------------- Plot (poster layout) ----------------------
col.OPL = [52 166 213]/255;
col.IPL = [83 206 144]/255;
col.GCL = [94 71 181]/255;
col.R = [252 193 97]/255;
col.B = [241 99 85]/255;

base.OPL = nanmean(OPL_data(1,:));
base.IPL = nanmean(IPL_data(1,:));
base.GCL = nanmean(GCL_data(1,:));

figure('Color','w','Position',[80 80 1100 600]);

% OPL
ax1 = subplot(2,3,1); hold on;
plot_band(tgrid, qO_lo, qO_hi, col.OPL, 0.15);
plot(tgrid, qO_md, 'Color', col.OPL, 'LineWidth', 2);
plot_layer_mean_sd(tdata, OPL_data, col.OPL, base.OPL);
title('Outer Plexiform Layer (OPL)','FontWeight','bold');
xlabel('Time (h)'); ylabel('Microglia density (cells/mm^2)'); fmt_axis();

% IPL
ax2 = subplot(2,3,2); hold on;
plot_band(tgrid, qI_lo, qI_hi, col.IPL, 0.15);
plot(tgrid, qI_md, 'Color', col.IPL, 'LineWidth', 2);
plot_layer_mean_sd(tdata, IPL_data, col.IPL, base.IPL);
title('Inner Plexiform Layer (IPL)','FontWeight','bold');
xlabel('Time (h)'); ylabel('Microglia density (cells/mm^2)'); fmt_axis();

% GCL/NFL
ax3 = subplot(2,3,3); hold on;
plot_band(tgrid, qG_lo, qG_hi, col.GCL, 0.15);
plot(tgrid, qG_md, 'Color', col.GCL, 'LineWidth', 2);
plot_layer_mean_sd(tdata, GCL_data, col.GCL, base.GCL);
title({'Ganglion Cell & Nerve Fiber Layers','(GCL/NFL)'},'FontWeight','bold');
xlabel('Time (h)'); ylabel('Microglia density (cells/mm^2)'); fmt_axis();

% ---- match y-limits AND tick positions across upper panels ----
yl = ylim(ax1);
yt = yticks(ax1);
ylim(ax2, yl); yticks(ax2, yt);
ylim(ax3, yl); yticks(ax3, yt);

% Peritoneum (R)
ax4 = subplot(2,3,4); hold on;
plot_band(tgrid, qR_lo, qR_hi, col.R, 0.12);
plot(tgrid, qR_md, 'Color', col.R, 'LineWidth', 2);
yline(0,':','Color',[.6 .6 .6],'LineWidth',1.2);
title('Peritoneum','FontWeight','bold');
xlabel('Time (h)'); ylabel('LPS Concentration (ng/ml)'); fmt_axis();

% Blood (B)
ax5 = subplot(2,3,5); hold on;
plot_band(tgrid, qB_lo, qB_hi, col.B, 0.12);
plot(tgrid, qB_md, 'Color', col.B, 'LineWidth', 2);
yline(0,':','Color',[.6 .6 .6],'LineWidth',1.2);
title('Blood','FontWeight','bold');
xlabel('Time (h)'); ylabel('LPS Concentration (ng/ml)'); fmt_axis();

% ---- match x-limits of upper panels to lower panels ----
xl = xlim(ax4);
xlim(ax1, xl); xlim(ax2, xl); xlim(ax3, xl);

% Blank
subplot(2,3,6); axis off;
sgtitle('','FontWeight','bold');

%% ---------------------- Helpers ----------------------
function p2 = overwrite_from_q(p2, q)
p2.eta_G = q(1); p2.k_OI = q(2); p2.k_IO = q(3); p2.k_IG = q(4); p2.k_G = q(5);
end

function [obsVec, timeVec, layerVec] = pack_observations(OPL_data, IPL_data, GCL_data, tdata)
obsVec=[]; timeVec=[]; layerVec=[];
    function append_layer(A, layer_id)
        for r = 1:numel(tdata)
            row = A(r,:); m = ~isnan(row);
            if any(m)
                obsVec   = [obsVec,   row(m)];
                timeVec  = [timeVec,  repmat(tdata(r),1,sum(m))];
                layerVec = [layerVec, repmat(layer_id,1,sum(m))];
            end
        end
    end
append_layer(OPL_data,1); append_layer(IPL_data,2); append_layer(GCL_data,3);
obsVec=obsVec(:); timeVec=timeVec(:); layerVec=layerVec(:);
end

function r = residuals_vec(q, obsVec, timeVec, layerVec, ICs, tspan, p)
p2 = overwrite_from_q(p, q);
tgrid   = linspace(tspan(1), tspan(2), 500);
opts_ode= odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',(tspan(2)-tspan(1))/10000);
try
    [t,y] = ode15s(@(tt,yy) microglia_model_opt12(tt,yy,p2), tgrid, ICs, opts_ode);
    m1 = interp1(t,y(:,1),timeVec,'linear','extrap');
    m2 = interp1(t,y(:,2),timeVec,'linear','extrap');
    m3 = interp1(t,y(:,3),timeVec,'linear','extrap');
    model = m1; model(layerVec==2)=m2(layerVec==2); model(layerVec==3)=m3(layerVec==3);
    bad = ~isfinite(model); model(bad)=0;
    r = obsVec - model; r(bad) = 1e6;
catch
    r = 1e6*ones(size(obsVec));
end
end

function lp = logpost_theta(theta, s, obsVec, timeVec, layerVec, ICs, tspan, p, LO, HI)
q = exp(theta(:)).';
sigma = exp(s);
if any(q < LO) || any(q > HI) || ~isfinite(sigma) || sigma<=0
    lp = -Inf; return;
end
r   = residuals_vec(q, obsVec, timeVec, layerVec, ICs, tspan, p);
SSE = sum(r.^2);
N   = numel(r);
lp  = -0.5*SSE/(sigma^2) - N*log(sigma) - 0.5*N*log(2*pi);
end

function plot_band(t, ylo, yhi, colorRGB, alphaVal)
t = t(:); ylo = ylo(:); yhi = yhi(:);
fill([t; flipud(t)], [ylo; flipud(yhi)], colorRGB, ...
    'FaceAlpha', alphaVal, 'EdgeColor', 'none');
end

function plot_layer_mean_sd(tdata, dataMat, colorRGB, baselineVal)
mn  = mean(dataMat, 2, 'omitnan');
sd  = std(dataMat,  0, 2, 'omitnan');
xr  = max(tdata) - min(tdata);
cap = max(0.02 * xr, 0.5);
lw  = 1.0;
ms  = 6;
plot(tdata, mn, 'o', 'MarkerSize', ms, 'LineWidth', lw, ...
    'MarkerEdgeColor', colorRGB, 'MarkerFaceColor', 'w', 'Color', colorRGB);
for i = 1:numel(tdata)
    line([tdata(i) tdata(i)], [mn(i)-sd(i), mn(i)+sd(i)], ...
        'Color', colorRGB, 'LineWidth', lw);
    line([tdata(i)-cap tdata(i)+cap], repmat(mn(i)-sd(i),1,2), ...
        'Color', colorRGB, 'LineWidth', lw);
    line([tdata(i)-cap tdata(i)+cap], repmat(mn(i)+sd(i),1,2), ...
        'Color', colorRGB, 'LineWidth', lw);
end
if isfinite(baselineVal)
    yline(baselineVal,'--','Color',[0.75 0.75 0.75],'LineWidth',1.2);
end
end

function fmt_axis()
set(gca,'LineWidth',2,'FontWeight','bold','TickDir','out'); box off;
end