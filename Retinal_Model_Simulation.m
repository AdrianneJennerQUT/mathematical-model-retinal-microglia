%% Retina_Model_Simulation.m
% Simulates and plots 6 panels in the screenshot order:
%   Top row:    (A) PK fit (Blood) | (B) Peritoneum | (C) Blood
%   Bottom row: (D) OPL            | (E) IPL        | (F) GCL/NFL
%
% Self-contained in ONE file (this script). Requires only Retinal_data.mat.

clear; clc;

%% ===================== Styling / Thickness Hierarchy =====================
LW.model = 3.5;   % model trajectories (dominant)
LW.mean  = 2.0;   % mean markers (edge)
LW.base  = 1.2;   % baseline / reference lines (subtle)
LW.rep   = 0.8;   % replicate markers (faint)
LW.axis  = 2.0;   % axes spines

MS.rep  = 4;      % replicate marker size
MS.mean = 8;      % mean marker size

%% ---------------------- Load data safely ----------------------
matfile = 'Retinal_data.mat';
if ~exist(matfile,'file')
    error('File "%s" not found in %s', matfile, pwd);
end
S = load(matfile);

required = {'time','OPL_individual','IPL_individual','GCL_NFL_individual'};
missing = required(~isfield(S, required));
if ~isempty(missing)
    error('Missing variable(s) in %s: %s', matfile, strjoin(missing, ', '));
end

tdata    = S.time(:);
OPL_data = S.OPL_individual;
IPL_data = S.IPL_individual;
GCL_data = S.GCL_NFL_individual;

rowmean  = @(A) mean(A,2,'omitnan');
OPL_mean = rowmean(OPL_data);
IPL_mean = rowmean(IPL_data);
GCL_mean = rowmean(GCL_data);

%% ---------------------- Microglia model parameters ----------------------
p.k_R   = 1.71;      % 1/h
p.F     = 0.03;      % dimensionless
p.Vol_P = 0.06;      % mL
p.Vol_B = 0.9548;    % mL
p.k_B   = 0.05776;   % 1/hm

% Half-activation (same units as B, ng/mL)
p.eta_O = 1667;
p.eta_G = 210;

% Fitted kinetics (1/h)
p.k_OI  = 0.221192;
p.k_IO  = 0.42835;
p.k_IG  = 0.0572969;
p.k_G   = 0.0610402;

%% ---------------------- Initial conditions & time grid ----------------------
% States: [O, I, G, R, B]
ICs   = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2.4889e+05*9, 0];

% Baseline G(0) for the logistic k_G term
p.G0 = ICs(3);

tspan = [tdata(1), tdata(end)];
tgrid = linspace(tspan(1), tspan(2), 800);

%% ---------------------- Simulate microglia ODE ----------------------
opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',(tspan(2)-tspan(1))/500);
[t_sim, y_sim] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p), tgrid, ICs, opts_ode);

O_sim = y_sim(:,1);
I_sim = y_sim(:,2);
G_sim = y_sim(:,3);
R_sim = y_sim(:,4);
B_sim = y_sim(:,5);

%% ---------------------- Colors ----------------------
col.OPL = [52 166 213]/255;
col.IPL = [83 206 144]/255;
col.GCL = [94 71 181]/255;
col.R   = [252 193 97]/255;
col.B   = [241 99 85]/255;

base.OPL = nanmean(OPL_data(1,:));
base.IPL = nanmean(IPL_data(1,:));
base.GCL = nanmean(GCL_data(1,:));

%% ---------------------- Figure layout (2x3) - screenshot order ----------------------
fig = figure('Color','w','Position',[80 80 1100 600]);
tlo = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');

% ===================== TOP ROW =====================

% (A) PK fit panel (Blood)
axA = nexttile(tlo, 1); hold(axA,'on');
plot_pk_fit_panel(axA, LW);

% (B) Peritoneum (R)
axB = nexttile(tlo, 2); hold(axB,'on');
plot(axB, t_sim, R_sim, 'Color', col.R, 'LineWidth', LW.model);
yline(axB, 0, ':', 'Color',[0.6 0.6 0.6], 'LineWidth', LW.base);
fmt_axis(axB, LW);
xlabel(axB,'Time (h)');
ylabel(axB,'LPS Concentration (ng/ml)');
title(axB,'Peritoneum','FontWeight','bold');

% (C) Blood (B)
axC = nexttile(tlo, 3); hold(axC,'on');
plot(axC, t_sim, B_sim, 'Color', col.B, 'LineWidth', LW.model);
yline(axC, 0, ':', 'Color',[0.6 0.6 0.6], 'LineWidth', LW.base);
fmt_axis(axC, LW);
xlabel(axC,'Time (h)');
ylabel(axC,'LPS Concentration (ng/ml)');
title(axC,'Blood','FontWeight','bold');

% ===================== BOTTOM ROW =====================

% (D) OPL
axD = nexttile(tlo, 4); hold(axD,'on');
plot_layer_with_means(axD, tdata, OPL_data, O_sim, t_sim, col.OPL, base.OPL, LW, MS);
title(axD,'Outer Plexiform Layer (OPL)','FontWeight','bold');

% (E) IPL
axE = nexttile(tlo, 5); hold(axE,'on');
plot_layer_with_means(axE, tdata, IPL_data, I_sim, t_sim, col.IPL, base.IPL, LW, MS);
title(axE,'Inner Plexiform Layer (IPL)','FontWeight','bold');

% (F) GCL/NFL
axF = nexttile(tlo, 6); hold(axF,'on');
plot_layer_with_means(axF, tdata, GCL_data, G_sim, t_sim, col.GCL, base.GCL, LW, MS);
title(axF,'Ganglion Cell Layer & Nerve Fiber Layer (GCL & NFL)','FontWeight','bold');

drawnow;

%% ====================== Local functions (must be below script) ======================

function dydt = rhs_microglia(~, y, p)
    O = y(1);
    I = y(2);
    G = y(3);
    R = y(4);
    B = y(5);

    hill_G = B / (B + p.eta_G);
    hill_O = B / (B + p.eta_O);

    % Logistic-style k_G * G * (1 - G/G0) flux between I and G
    flux_G = p.k_G * G * (1 - G / p.G0);

    dO = -p.k_OI * O * hill_G + p.k_IO * I * hill_O;

    dI =  p.k_OI * O * hill_G ...
        - p.k_IO * I * hill_O ...
        - p.k_IG * I * hill_G ...
        - flux_G;

    dG =  p.k_IG * I * hill_G ...
        + flux_G;

    dR = -p.k_R * R;

    dB =  p.F * p.k_R * (p.Vol_P / p.Vol_B) * R ...
        - p.k_B * B;

    dydt = [dO; dI; dG; dR; dB];
end

function plot_layer_with_means(ax, tdata, dataMat, y_sim, t_sim, colorRGB, baselineVal, LW, MS)
    % Replicate points (small open circles, faint)
    light = colorRGB + (1 - colorRGB)*0.60;
    [T, ~] = size(dataMat);
    for r = 1:T
        row = dataMat(r,:); m = ~isnan(row);
        if any(m)
            plot(ax, tdata(r)*ones(1,sum(m)), row(m), 'o', ...
                'MarkerSize', MS.rep, ...
                'Color', light, ...
                'LineWidth', LW.rep, ...
                'MarkerFaceColor','none');
        end
    end

    % Mean at each time (large white circle with colored edge)
    means = mean(dataMat,2,'omitnan');
    plot(ax, tdata, means, 'o', ...
        'MarkerSize', MS.mean, ...
        'LineWidth', LW.mean, ...
        'MarkerFaceColor','w', ...
        'MarkerEdgeColor', colorRGB);

    % Model curve (dominant)
    plot(ax, t_sim, y_sim, 'Color', colorRGB, 'LineWidth', LW.model);

    % Baseline (dashed grey, subtle)
    if isfinite(baselineVal)
        yline(ax, baselineVal,'--','Color',[0.75 0.75 0.75],'LineWidth', LW.base);
    end

    fmt_axis(ax, LW);
    xlabel(ax,'Time (h)');
    ylabel(ax,'Microglia density (cells/mm^2)');
end

function fmt_axis(ax, LW)
    set(ax,'LineWidth',LW.axis,'FontWeight','bold','TickDir','out');
    box(ax,'off');
    try
        ax.Toolbar.Visible = 'off';
    catch
    end
end

%% ====================== PK fitting panel (Python colors + legend) ======================

function plot_pk_fit_panel(ax, LW)
    % Data points
    time_points = [0.5, 2, 4, 12, 24]';
    dose_0_5 = [136.986, 223.287, 206.849, 98.630, 60.274]';
    dose_1_0 = [250.685, 372.603, 363.014, 242.466, 139.726]';
    combined = [dose_0_5; dose_1_0];

    % Error bars (0.5 mg/kg)
    high_05 = [169.15, 260.64, 222.34, 119.15, 78.72]';
    low_05  = [108.51, 192.55, 190.43, 81.91, 41.49]';
    err05_lower = dose_0_5 - low_05;
    err05_upper = high_05 - dose_0_5;

    % Error bars (1.0 mg/kg)
    high_10 = [274.47, 398.94, 394.68, 282.98, 152.13]';
    low_10  = [219.15, 350, 345.74, 202.13, 125.53]';
    err10_lower = dose_1_0 - low_10;
    err10_upper = high_10 - dose_1_0;

    % Fixed params (match Python)
    Vol_P = 0.06;
    Vol_B = 0.9548;
    k2    = 0.05776;

    % Exact Python colors
    col_data_05  = [1, 0.5, 0.5];
    col_data_10  = [0.5, 0.2, 0.9];
    col_model_05 = [1, 0.75, 0.8];
    col_model_10 = [0.7, 0.4, 0.9];
    ci_alpha = 0.30;

    % Fit in log-space (matches Python)
    params0 = [1, 0.05];  % [kPB, F]
    obj = @(p) (FObjective_pk(p, time_points, Vol_P, Vol_B, k2) - log(combined));

    use_lsqnonlin = exist('lsqnonlin','file') == 2;
    if use_lsqnonlin
        opts = optimoptions('lsqnonlin', ...
            'Display','off', ...
            'FiniteDifferenceType','central', ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12);
        [popt, ~, resid, ~, ~, ~, J] = lsqnonlin(obj, params0, [], [], opts);
    else
        popt  = fminsearch(@(p) sum(obj(p).^2), params0);
        resid = obj(popt);
        J     = finiteDiffJacobian_pk(obj, popt);
    end

    kPB_hat = popt(1);
    F_hat   = popt(2);

    % Covariance estimate -> 95% CI like Python (±1.96*SE)
    N   = numel(resid);
    np  = numel(popt);
    SSE = sum(resid.^2);
    s2  = SSE / max(N - np, 1);
    pcov = s2 * inv(J.'*J); %#ok<MINV>
    perr = sqrt(diag(pcov));

    z = 1.96;
    kPB_lo = kPB_hat - z*perr(1);
    kPB_hi = kPB_hat + z*perr(1);
    F_lo   = F_hat   - z*perr(2);
    F_hi   = F_hat   + z*perr(2);

    % Curves
    tfine = linspace(0, 24, 1000)';

    model_05 = getModelOutput_pk([kPB_hat, F_hat], tfine, 0.5, Vol_P, Vol_B, k2);
    model_10 = getModelOutput_pk([kPB_hat, F_hat], tfine, 1.0, Vol_P, Vol_B, k2);

    lo_05 = getModelOutput_pk([kPB_lo, F_lo], tfine, 0.5, Vol_P, Vol_B, k2);
    hi_05 = getModelOutput_pk([kPB_hi, F_hi], tfine, 0.5, Vol_P, Vol_B, k2);

    lo_10 = getModelOutput_pk([kPB_lo, F_lo], tfine, 1.0, Vol_P, Vol_B, k2);
    hi_10 = getModelOutput_pk([kPB_hi, F_hi], tfine, 1.0, Vol_P, Vol_B, k2);

    % Plot CI bands first (legend handles)
    hCI05 = fill_ci(ax, tfine, lo_05, hi_05, col_model_05, ci_alpha);
    hCI10 = fill_ci(ax, tfine, lo_10, hi_10, col_model_10, ci_alpha);

    % Plot model lines
    hM05  = plot(ax, tfine, model_05, '-', 'Color', col_model_05, 'LineWidth', 1.8);
    hM10  = plot(ax, tfine, model_10, '-', 'Color', col_model_10, 'LineWidth', 1.8);

    % Plot data + error bars (marker handles for legend)
    hD05  = plot_asym_err(ax, time_points, dose_0_5, err05_lower, err05_upper, col_data_05, 6);
    hD10  = plot_asym_err(ax, time_points, dose_1_0, err10_lower, err10_upper, col_data_10, 6);

    xlabel(ax,'Time (h)');
    ylabel(ax,'LPS Concentration (ng/mL)');
    title(ax,'PK Subsystem Fitting (Blood)','FontWeight','bold');

    grid(ax,'off');
    xticks(ax, 0:2:24);
    xlim(ax, [-1 25]);

    fmt_axis(ax, LW);

    % Legend matching your screenshot order
    leg = legend(ax, [hD05, hD10, hM05, hM10, hCI05, hCI10], ...
        {'Observed data 0.5 mg/kg', ...
         'Observed data 1.0 mg/kg', ...
         'Model (0.5 mg/kg)', ...
         'Model (1.0 mg/kg)', ...
         'Model 95% CI (0.5 mg/kg)', ...
         'Model 95% CI (1.0 mg/kg)'}, ...
        'Location','northeast', 'Box','on');
    leg.FontSize = 6.5;
end

function ylog = FObjective_pk(params, time_points, Vol_P, Vol_B, k2)
    B05 = getModelOutput_pk(params, time_points, 0.5, Vol_P, Vol_B, k2);
    B10 = getModelOutput_pk(params, time_points, 1.0, Vol_P, Vol_B, k2);
    epsB = 1e-12;
    ylog = [log(max(B05, epsB)); log(max(B10, epsB))];
end

function B_sol = getModelOutput_pk(params, t_eval, dose, Vol_P, Vol_B, k2)
    kPB = params(1);
    F   = params(2);

    P0 = dose * 248000;
    B0 = 0;
    y0 = [P0; B0];

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t, y] = ode45(@(t,y) systemODE_pk(t,y,kPB,F,Vol_P,Vol_B,k2), [0 24], y0, opts);

    B_sol = interp1(t, y(:,2), t_eval, 'pchip');
end

function dydt = systemODE_pk(~, Y, kPB, F, Vol_P, Vol_B, k2)
    P = Y(1);
    B = Y(2);
    dPdt = -kPB * P;
    dBdt = F * kPB * (Vol_P / Vol_B) * P - k2 * B;
    dydt = [dPdt; dBdt];
end

function J = finiteDiffJacobian_pk(fun, p0)
    m  = numel(fun(p0));
    n  = numel(p0);
    J  = zeros(m,n);

    h = 1e-6 * max(1, abs(p0(:)));
    for j = 1:n
        p1 = p0; p2 = p0;
        p1(j) = p1(j) + h(j);
        p2(j) = p2(j) - h(j);
        J(:,j) = (fun(p1) - fun(p2)) / (2*h(j));
    end
end

function h = fill_ci(ax, x, ylo, yhi, col, alphaVal)
    x = x(:); ylo = ylo(:); yhi = yhi(:);
    h = fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], col, ...
        'LineStyle','none', 'EdgeColor','none');
    set(h,'FaceAlpha', alphaVal);
    uistack(h,'bottom');
end

function h = plot_asym_err(ax, x, y, elow, ehigh, col, markersize)
    xr  = max(x) - min(x);
    cap = max(0.02 * xr, 0.15);

    h = plot(ax, x, y, 'o', ...
        'MarkerSize', markersize, ...
        'LineWidth', 1.0, ...
        'MarkerEdgeColor', col, ...
        'MarkerFaceColor', 'w', ...
        'Color', col);

    for i = 1:numel(x)
        line(ax, [x(i) x(i)], [y(i)-elow(i), y(i)+ehigh(i)], 'Color', col, 'LineWidth', 1.0);
        line(ax, [x(i)-cap x(i)+cap], [y(i)-elow(i) y(i)-elow(i)], 'Color', col, 'LineWidth', 1.0);
        line(ax, [x(i)-cap x(i)+cap], [y(i)+ehigh(i) y(i)+ehigh(i)], 'Color', col, 'LineWidth', 1.0);
    end
end