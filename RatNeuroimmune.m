function Retinal_Model_3x3_Combined
clear; clc; close all;

%% ---------------------- Load data safely ----------------------
matfile = 'Retinal_data.mat';
if ~exist(matfile,'file')
    error('File "%s" not found in %s', matfile, pwd);
end
S = load(matfile);

required = {'time','OPL_individual','IPL_individual','GCL_NFL_individual'};
missing  = required(~isfield(S, required));
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

%% ---------------------- Colors ----------------------
col.OPL = [52  166 213]/255;
col.IPL = [83  206 144]/255;
col.GCL = [94  71  181]/255;
col.R   = [252 193  97]/255;
col.B   = [241  99  85]/255;

%% =====================================================================
%  PARAMETERS
%% =====================================================================

p_mouse.k_R   = 1.71;
p_mouse.F     = 0.03;
p_mouse.Vol_P = 0.06;
p_mouse.Vol_B = 0.9548;
p_mouse.k_B   = 0.05776;
p_mouse.eta_O = 1667;
p_mouse.eta_G = 210;
p_mouse.k_OI  = 0.221192;
p_mouse.k_IO  = 0.42835;
p_mouse.k_IG  = 0.0572969;
p_mouse.k_G   = 0.0610402;
p_mouse.G0    = GCL_mean(1);

p_rat.k_R   = 0.126442;
p_rat.F     = 0.0013;
p_rat.Vol_P = 3.07;
p_rat.Vol_B = 18.112;
p_rat.k_B   = 0.0933;
p_rat.eta_O = 1667;
p_rat.eta_G = 210;
p_rat.k_OI  = 0.221192;
p_rat.k_IO  = 0.42835;
p_rat.k_IG  = 0.0572969;
p_rat.k_G   = 0.0610402;
p_rat.G0    = GCL_mean(1);

%% ---------------------- Time grids ----------------------
tspan    = [tdata(1), tdata(end)];
tgrid    = linspace(tspan(1), tspan(2), 800);
opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',(tspan(2)-tspan(1))/500);

%% =====================================================================
%  ROW 1 PANEL A: PK Subsystem Fitting
%% =====================================================================

t_fixed_g = [0.25 0.5 1 2 4 6 8 12 18 24 48]';
c_GN_g = [ ...
    0.08974359; 0.192307692; 0.397435897; ...
    1.474358974358974; 1.2435897435897434; 1.782051282051282; ...
    1.9871794871794872; 1.801282051282051; 1.3141025641025639; ...
    0.4294871794871789; 0.17948717948717746 ...
] * 1000;

sq_data_g = [ ...
    0.044871795 0.134615385;
    0.173076923 0.294871795;
    0.58974359  0.294871795;
    1.788461538 1.16025641;
    1.378205128 1.108974359;
    1.987179487 1.583333333;
    2.820512821 1.153846154;
    1.993589744 1.602564103;
    1.455128205 1.179487179;
    0.544871795 0.294871795;
    0.147435897 0.211538462 ...
] * 1000;

Vol_P_g = 3.07;
Vol_B_g = 18.112;
k2_g    = 0.0933;
dose_g  = 0.2;
P0_g    = dose_g * 92000;
B0_g    = 0;

x0_g  = [1.0, 0.05];
lb_g  = [0, 0];
ub_g  = [Inf, Inf];

optsFit = optimoptions('lsqcurvefit','Display','off', ...
    'FiniteDifferenceType','central','FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12,'MaxIterations',2000,'MaxFunctionEvaluations',20000);

modelFun_g = @(x,t) model_B_ng_per_mL(t, x(1), x(2), Vol_P_g, Vol_B_g, k2_g, P0_g, B0_g);

log_c_GN_g = log(c_GN_g);

modelFun_g_log = @(x,t) log( ...
    max(model_B_ng_per_mL(t, x(1), x(2), ...
    Vol_P_g, Vol_B_g, k2_g, P0_g, B0_g), eps) );

[xhat_g, ~, resid_g, ~, ~, ~, J_g] = lsqcurvefit( ...
    modelFun_g_log, x0_g, t_fixed_g, log_c_GN_g, lb_g, ub_g, optsFit);

kPB_hat = xhat_g(1);
F_hat   = xhat_g(2);

N_g   = numel(resid_g);
p_g   = numel(xhat_g);
SSE_g = sum(resid_g.^2);
s2_g  = SSE_g / max(N_g - p_g, 1);
pcov_g = s2_g * inv(J_g.'*J_g); %#ok<MINV>
perr_g = sqrt(diag(pcov_g));

z        = 1.96;
kPB_low  = kPB_hat - z*perr_g(1);
kPB_high = kPB_hat + z*perr_g(1);
F_low    = F_hat   - z*perr_g(2);
F_high   = F_hat   + z*perr_g(2);

y_pred_g = modelFun_g(xhat_g, t_fixed_g);
rmse_g   = sqrt(mean((c_GN_g - y_pred_g).^2));

%% ---- Print Panel A parameters ----
fprintf('\n============================================\n');
fprintf('PANEL A: PK SUBSYSTEM FIT PARAMETERS\n');
fprintf('============================================\n');

fprintf('Fitted parameters:\n');
fprintf('  kPB_hat   = %.6f\n', kPB_hat);
fprintf('  F_hat     = %.6f\n', F_hat);

fprintf('\n95%% confidence interval bounds:\n');
fprintf('  kPB_low   = %.6f\n', kPB_low);
fprintf('  kPB_high  = %.6f\n', kPB_high);
fprintf('  F_low     = %.6f\n', F_low);
fprintf('  F_high    = %.6f\n', F_high);

fprintf('\nFixed PK parameters used:\n');
fprintf('  Vol_P_g   = %.6f mL\n', Vol_P_g);
fprintf('  Vol_B_g   = %.6f mL\n', Vol_B_g);
fprintf('  k2_g      = %.6f 1/h\n', k2_g);
fprintf('  dose_g    = %.6f mg/kg\n', dose_g);
fprintf('  P0_g      = %.6f\n', P0_g);
fprintf('  B0_g      = %.6f\n', B0_g);

fprintf('\nLeast-squares fitting setup:\n');
fprintf('  Initial guess x0_g = [%.6f, %.6f]\n', x0_g(1), x0_g(2));
fprintf('  Lower bounds lb_g  = [%.6f, %.6f]\n', lb_g(1), lb_g(2));
fprintf('  Upper bounds ub_g  = [%.6f, %.6f]\n', ub_g(1), ub_g(2));

fprintf('\nFit diagnostics:\n');
fprintf('  SSE       = %.6f\n', SSE_g);
fprintf('  RMSE      = %.6f ng/mL\n', rmse_g);
fprintf('  N points  = %d\n', N_g);

fprintf('\nObserved data used for fitting:\n');
disp(table(t_fixed_g, c_GN_g, ...
    'VariableNames', {'Time_h','Observed_LPS_ng_per_mL'}));

fprintf('============================================\n\n');

% Extended from 20 h to 50 h
t_plot_g = linspace(0, 50, 1000)';
B_best_g = model_B_ng_per_mL(t_plot_g, kPB_hat, F_hat,   Vol_P_g, Vol_B_g, k2_g, P0_g, B0_g);
B_low_g  = model_B_ng_per_mL(t_plot_g, kPB_low,  F_low,  Vol_P_g, Vol_B_g, k2_g, P0_g, B0_g);
B_high_g = model_B_ng_per_mL(t_plot_g, kPB_high, F_high, Vol_P_g, Vol_B_g, k2_g, P0_g, B0_g);

%% =====================================================================
%  ROW 2: Mouse vs Rat model simulations
%% =====================================================================

IC_mouse = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2.4889e+05*9, 0];
IC_rat   = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 828000, 0];

[t_mouse, y_mouse] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p_mouse), tgrid, IC_mouse, opts_ode);
[t_rat,   y_rat  ] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p_rat),   tgrid, IC_rat,   opts_ode);

O_mouse = y_mouse(:,1); I_mouse = y_mouse(:,2); G_mouse = y_mouse(:,3);
R_mouse = y_mouse(:,4); B_mouse = y_mouse(:,5);

O_rat   = y_rat(:,1);   I_rat   = y_rat(:,2);   G_rat   = y_rat(:,3);
R_rat   = y_rat(:,4);   B_rat   = y_rat(:,5);

%% =====================================================================
%  ROW 3: Dose sweep
%% =====================================================================

doseMults = [0 0.25 0.5 1 2 3 5 7 10];
nMult     = numel(doseMults);

cmap_full = readmatrix('deep.txt');
cmap = interp1(linspace(0,1,size(cmap_full,1)), cmap_full, linspace(0,1,nMult));
cmap(1,:) = [0.99, 0.94, 0.70];

%% =====================================================================
%  FIGURE LAYOUT
%% =====================================================================

LW      = 2.0;
LW_BASE = 1.5;
titleFS = 10.5;

left  = 0.07;
right = 0.03;
hgap  = 0.06;
vgap  = 0.105;

axW = (1 - left - right - 2*hgap) / 3;
axH = 0.22;

row1_bot = 0.72;
row2_bot = row1_bot - axH - vgap;
row3_bot = row2_bot - axH - vgap;

x1 = left;
x2 = left + axW + hgap;
x3 = left + 2*(axW + hgap);

fig = figure('Color','w','Position',[60 60 900 900]);

%% =====================================================================
%  ROW 1 AXES
%% =====================================================================

axA = axes('Parent',fig,'Position',[x1 row1_bot axW axH]); hold(axA,'on');
axB = axes('Parent',fig,'Position',[x2 row1_bot axW axH]); hold(axB,'on');
axC = axes('Parent',fig,'Position',[x3 row1_bot axW axH]); hold(axC,'on');

%% ---- Panel A: PK fit ----

col_model_g = [0.2 0.8 0.8];
col_err_g   = [93 58 155]/255;
ci_alpha_g  = 0.30;

hCI = fill(axA, [t_plot_g; flipud(t_plot_g)], [B_low_g; flipud(B_high_g)], ...
    col_model_g, 'LineStyle','none', 'FaceAlpha', ci_alpha_g);
uistack(hCI,'bottom');

hM = plot(axA, t_plot_g, B_best_g, '-', 'Color', col_model_g, 'LineWidth', 1.8);

% Include all data points up to 50 h, including the 48 h point
mask = t_fixed_g <= 50;
t_plot_pts  = t_fixed_g(mask);
c_plot_pts  = c_GN_g(mask);
sq_plot_pts = sq_data_g(mask,:);

hD = plot(axA, t_plot_pts, c_plot_pts, 'o', 'MarkerSize',6, 'LineWidth',1.0, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineStyle','none');

capw = 0.5;
for i = 1:numel(t_plot_pts)
    t  = t_plot_pts(i);
    lo = sq_plot_pts(i,1);
    hi = sq_plot_pts(i,2);
    line(axA, [t t],          [lo hi], 'Color', col_err_g, 'LineWidth', 1.4);
    line(axA, [t-capw t+capw],[hi hi], 'Color', col_err_g, 'LineWidth', 1.4);
    line(axA, [t-capw t+capw],[lo lo], 'Color', col_err_g, 'LineWidth', 1.4);
end

xlabel(axA,'Time (hours)');
ylabel(axA,'LPS Concentration (ng/mL)');
title(axA,'PK Subsystem Fitting (Blood)','FontWeight','bold','FontSize',titleFS);

xlim(axA,[0 50]);
xticks(axA,0:10:50);
ylim(axA,[0 3000]);

% Legend-style box for data/model/CI
hObs_leg = plot(axA, NaN, NaN, 'o', ...
    'MarkerSize',6, 'LineWidth',1.2, ...
    'MarkerEdgeColor','k', 'MarkerFaceColor','w');

hCurve_leg = plot(axA, NaN, NaN, '-', ...
    'Color', col_model_g, 'LineWidth',2.0);

hShade_leg = patch(axA, NaN, NaN, col_model_g, ...
    'FaceAlpha',ci_alpha_g, 'EdgeColor','none');

legA = legend(axA, [hObs_leg hCurve_leg hShade_leg], ...
    {'Observed data 0.2 mg/kg', ...
     'Model (0.2 mg/kg)', ...
     'Model 95% CI (0.2 mg/kg)'}, ...
    'Location','none', ...
    'Box','on', ...
    'FontSize',7.5, ...
    'FontWeight','bold');

legA.LineWidth = 1.5;

% [left bottom width height]
legA.Position = [0.135 0.88 0.19 0.050];

fmt_axis(axA);

% RMSE-only box underneath the legend
txtA = sprintf('RMSE = %.3f', rmse_g);

annotation(fig,'textbox', [x1 + 0.202, row1_bot + 0.150, 0.001, 0.001], ...
    'String', txtA, ...
    'FitBoxToText','on', ...
    'Margin', 1.8, ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'LineWidth',1.5, ...
    'FontWeight','bold', ...
    'FontSize',8, ...
    'HorizontalAlignment','center');
%% ---- Panel B: Peritoneum ----

plot(axB, t_mouse, R_mouse, '-',  'Color', col.R, 'LineWidth', LW);
plot(axB, t_rat,   R_rat,   '--', 'Color', col.R, 'LineWidth', LW);
xlabel(axB,'Time (h)');
ylabel(axB,'LPS Concentration (ng/mL)');
title(axB,'Peritoneum','FontWeight','bold','FontSize',titleFS);
fmt_axis(axB);
legend(axB,{'Mouse PK model','Rat PK model'},'Location','northeast','Box','off','FontSize',8);

%% ---- Panel C: Blood ----

plot(axC, t_mouse, B_mouse*1000, '-',  'Color', col.B, 'LineWidth', LW);
plot(axC, t_rat,   B_rat*1000,   '--', 'Color', col.B, 'LineWidth', LW);
xlabel(axC,'Time (h)');
ylabel(axC,'LPS Concentration (ng/mL)');
title(axC,'Blood','FontWeight','bold','FontSize',titleFS);
fmt_axis(axC);
legend(axC,{'Mouse PK model','Rat PK model'},'Location','northeast','Box','off','FontSize',8);

%% =====================================================================
%  ROW 2 AXES
%% =====================================================================

axD = axes('Parent',fig,'Position',[x1 row2_bot axW axH]); hold(axD,'on');
axE = axes('Parent',fig,'Position',[x2 row2_bot axW axH]); hold(axE,'on');
axF = axes('Parent',fig,'Position',[x3 row2_bot axW axH]); hold(axF,'on');

%% ---- Panel D: OPL ----

plot(axD, t_mouse, O_mouse, '-',  'Color', col.OPL, 'LineWidth', LW);
plot(axD, t_rat,   O_rat,   '--', 'Color', col.OPL, 'LineWidth', LW);
yline(axD, OPL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);
xlabel(axD,'Time (h)'); ylabel(axD,'Microglia density (cells/mm^2)');
title(axD,'Outer Plexiform Layer (OPL)','FontWeight','bold','FontSize',titleFS);
fmt_axis(axD);
ylim(axD,[80, 184]);
legend(axD,{'Mouse PK model','Rat PK model'},'Location','northeast','Box','off','FontSize',8);

%% ---- Panel E: IPL ----

plot(axE, t_mouse, I_mouse, '-',  'Color', col.IPL, 'LineWidth', LW);
plot(axE, t_rat,   I_rat,   '--', 'Color', col.IPL, 'LineWidth', LW);
yline(axE, IPL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);
xlabel(axE,'Time (h)'); ylabel(axE,'Microglia density (cells/mm^2)');
title(axE,'Inner Plexiform Layer (IPL)','FontWeight','bold','FontSize',titleFS);
fmt_axis(axE);
ylim(axE,[100, 273]);
legend(axE,{'Mouse PK model','Rat PK model'},'Location','northeast','Box','off','FontSize',8);

%% ---- Panel F: GCL/NFL ----

plot(axF, t_mouse, G_mouse, '-',  'Color', col.GCL, 'LineWidth', LW);
plot(axF, t_rat,   G_rat,   '--', 'Color', col.GCL, 'LineWidth', LW);
yline(axF, GCL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);
xlabel(axF,'Time (h)'); ylabel(axF,'Microglia density (cells/mm^2)');
title(axF,{'Ganglion Cell Layer &','Nerve Fiber Layer (GCL & NFL)'},'FontWeight','bold','FontSize',titleFS);
fmt_axis(axF);
ylim(axF,[80, 180]);
legend(axF,{'Mouse PK model','Rat PK model'},'Location','northeast','Box','off','FontSize',8);

%% =====================================================================
%  ROW 3 AXES: Dose sweep
%% =====================================================================

axG  = axes('Parent',fig,'Position',[x1 row3_bot axW axH]); hold(axG,'on');
axH2 = axes('Parent',fig,'Position',[x2 row3_bot axW axH]); hold(axH2,'on');
axI  = axes('Parent',fig,'Position',[x3 row3_bot axW axH]); hold(axI,'on');

yline(axG,  OPL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);
yline(axH2, IPL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);
yline(axI,  GCL_mean(1), '--', 'Color',[.75 .75 .75], 'LineWidth', LW_BASE);

lgdEntries  = cell(1, nMult);
hdlLegSweep = gobjects(1, nMult);

for k = 1:nMult
    mult = doseMults(k);
    ICs  = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2.4889e+05 * mult, 0];

    [t_sim, y_sim] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p_mouse), tgrid, ICs, opts_ode);

    plot(axG,  t_sim, y_sim(:,1), 'Color', cmap(k,:), 'LineWidth', 1.8);
    plot(axH2, t_sim, y_sim(:,2), 'Color', cmap(k,:), 'LineWidth', 1.8);
    plot(axI,  t_sim, y_sim(:,3), 'Color', cmap(k,:), 'LineWidth', 1.8);

    lgdEntries{k}  = sprintf('%.2f mg/kg LPS', mult);
    hdlLegSweep(k) = plot(axI, NaN, NaN, 'Color', cmap(k,:), 'LineWidth', 2.0);
end

xlabel(axG,'Time (h)'); ylabel(axG,'Microglia density (cells/mm^2)');
title(axG,'Outer Plexiform Layer (OPL)','FontWeight','bold','FontSize',titleFS);
fmt_axis(axG);

xlabel(axH2,'Time (h)'); ylabel(axH2,'Microglia density (cells/mm^2)');
title(axH2,'Inner Plexiform Layer (IPL)','FontWeight','bold','FontSize',titleFS);
fmt_axis(axH2);

xlabel(axI,'Time (h)'); ylabel(axI,'Microglia density (cells/mm^2)');
title(axI,{'Ganglion Cell Layer &','Nerve Fiber Layer (GCL & NFL)'},'FontWeight','bold','FontSize',titleFS);
fmt_axis(axI);

%% ---- Horizontal legend below row 3 ----

leg_fontsize    = 7.19;
leg_shift_right = 0.045;
leg_shift_down  = 0.02;

drawnow;
axLeg = axes('Parent',fig,'Position',[left, 0.005, 1-left-right, 0.04],'Visible','off');
hold(axLeg,'on');

hdlDummy = gobjects(1,nMult);
for k = 1:nMult
    hdlDummy(k) = plot(axLeg, NaN, NaN, 'Color', cmap(k,:), 'LineWidth', 2.0);
end

legSweep = legend(axLeg, hdlDummy, lgdEntries, ...
    'Orientation','horizontal','NumColumns',nMult, ...
    'Location','north','Box','off','FontSize',leg_fontsize,'FontWeight','bold');

legSweep.AutoUpdate = 'off';
try
    legSweep.ItemTokenSize = [13 9];
catch
end

drawnow;
legPos = legSweep.Position;
legSweep.Position = [legPos(1) + leg_shift_right, ...
                     legPos(2) - leg_shift_down, ...
                     legPos(3), legPos(4)];

drawnow;

end

%% =====================================================================
%  Helper functions
%% =====================================================================

function dydt = rhs_microglia(~, y, p)
    O = y(1); I = y(2); G = y(3); R = y(4); B = y(5);

    hill_G = B / (B + p.eta_G);
    hill_O = B / (B + p.eta_O);
    flux_G = p.k_G * G * (1 - G / p.G0);

    dO =  -p.k_OI * O * hill_G + p.k_IO * I * hill_O;
    dI =   p.k_OI * O * hill_G - p.k_IO * I * hill_O - p.k_IG * I * hill_G - flux_G;
    dG =   p.k_IG * I * hill_G + flux_G;
    dR =  -p.k_R * R;
    dB =   p.F * p.k_R * (p.Vol_P / p.Vol_B) * R - p.k_B * B;

    dydt = [dO; dI; dG; dR; dB];
end

function B_ng = model_B_ng_per_mL(t_eval, kPB, F, Vol_P, Vol_B, k2, P0, B0)
    t_eval = t_eval(:);
    tspan  = [0, max(t_eval)];

    odefun = @(t,Y) [-kPB*Y(1); F*kPB*(Vol_P/Vol_B)*Y(1) - k2*Y(2)];

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    sol  = ode45(odefun, tspan, [P0; B0], opts);
    Y    = deval(sol, t_eval);

    B_ng = Y(2,:).' * 1000;
end

function fmt_axis(ax)
    set(ax,'LineWidth',1.5,'FontWeight','bold','TickDir','out');
    box(ax,'off');
    try
        ax.Toolbar.Visible = 'off';
    catch
    end
end