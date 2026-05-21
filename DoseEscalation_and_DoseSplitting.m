%% Retinal_Model_Simulation_6Panel_Combined.m
clear; clc;

%% ---------------------- User controls ----------------------
titleFS = 11.5;

left   = 0.07;
right  = 0.06;
hgap   = 0.05;
axW    = (1 - left - right - 2*hgap) / 3;
topY   = 0.59;
botY   = 0.10;
axH    = 0.34;

legX   = 0.09;
legW   = 1 - left - right;
legY   = 0.48;
legH   = 0.055;

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

%% ---------------------- Parameters ----------------------
p.k_R   = 1.71;
p.F     = 0.03;
p.Vol_P = 0.06;
p.Vol_B = 0.9548;
p.k_B   = 0.0578;
p.eta_O = 1667;
p.eta_G = 210;
p.k_OI  = 0.221192;
p.k_IO  = 0.42835;
p.k_IG  = 0.0572969;
p.k_G   = 0.0610402;
p.G0    = GCL_mean(1);

%% ---------------------- Baselines ----------------------
base.OPL = nanmean(OPL_data(1,:));
base.IPL = nanmean(IPL_data(1,:));
base.GCL = nanmean(GCL_data(1,:));

%% ---------------------- Time grid ----------------------
tspan    = [tdata(1), tdata(end)];
tgrid    = linspace(tspan(1), tspan(2), 800);
opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',(tspan(2)-tspan(1))/500);

%% ---------------------- Dose sweep colormap ----------------------
doseMults = [0 0.25 0.5 1 2 3 5 7 10];
nMult     = numel(doseMults);
cmap_full = readmatrix('deep.txt');
cmap = interp1(linspace(0,1,size(cmap_full,1)), cmap_full, linspace(0,1,nMult));
cmap(1,:) = [1.00, 0.95, 0.60];

%% ---------------------- Split-dose simulations ----------------------
col.OPL = [52 166 213]/255;
col.IPL = [83 206 144]/255;
col.GCL = [94 71 181]/255;

IC_single = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2021000, 0];
R0_total  = IC_single(4);

[t_single, y_single] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p), tgrid, IC_single, opts_ode);
O_single = y_single(:,1);
I_single = y_single(:,2);
G_single = y_single(:,3);

[t_6,  y_6]  = simulate_split_peritoneum(p, IC_single, R0_total, tspan, 6,  opts_ode);
[t_12, y_12] = simulate_split_peritoneum(p, IC_single, R0_total, tspan, 12, opts_ode);
[t_24, y_24] = simulate_split_peritoneum(p, IC_single, R0_total, tspan, 24, opts_ode);
O_6  = y_6(:,1);  I_6  = y_6(:,2);  G_6  = y_6(:,3);
O_12 = y_12(:,1); I_12 = y_12(:,2); G_12 = y_12(:,3);
O_24 = y_24(:,1); I_24 = y_24(:,2); G_24 = y_24(:,3);

%% ---------------------- Figure layout ----------------------
fig = figure('Color','w','Position',[60 80 1200 820]);
x1 = left;
x2 = left + axW + hgap;
x3 = left + 2*(axW + hgap);

%% ---------------------- Top row: dose sweep ----------------------
ax1 = axes('Parent',fig,'Position',[x1 topY axW axH]); hold(ax1,'on');
fmt_axis(ax1); xlabel(ax1,'Time (h)'); ylabel(ax1,'Microglia density (cells/mm^2)');
title(ax1,'Outer Plexiform Layer (OPL)','FontWeight','bold','FontSize',titleFS);

ax2 = axes('Parent',fig,'Position',[x2 topY axW axH]); hold(ax2,'on');
fmt_axis(ax2); xlabel(ax2,'Time (h)'); ylabel(ax2,'Microglia density (cells/mm^2)');
title(ax2,'Inner Plexiform Layer (IPL)','FontWeight','bold','FontSize',titleFS);

ax3 = axes('Parent',fig,'Position',[x3 topY axW axH]); hold(ax3,'on');
fmt_axis(ax3); xlabel(ax3,'Time (h)'); ylabel(ax3,'Microglia density (cells/mm^2)');
title(ax3,'Ganglion Cell Layer & Nerve Fiber Layer (GCL & NFL)','FontWeight','bold','FontSize',titleFS);

lgdEntriesTop = cell(1, nMult);
hdlLegTop     = gobjects(1, nMult);

for k = 1:nMult
    mult = doseMults(k);
    ICs  = [OPL_mean(1), IPL_mean(1), GCL_mean(1), 2021000 * mult, 0];
    [t_sim, y_sim] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p), tgrid, ICs, opts_ode);
    plot(ax1, t_sim, y_sim(:,1), 'Color', cmap(k,:), 'LineWidth', 1.8);
    plot(ax2, t_sim, y_sim(:,2), 'Color', cmap(k,:), 'LineWidth', 1.8);
    plot(ax3, t_sim, y_sim(:,3), 'Color', cmap(k,:), 'LineWidth', 1.8);
    lgdEntriesTop{k} = sprintf('%.2f mg/kg LPS', mult);
end

% ---- Baseline dashed lines for top row ----
yline(ax1, base.OPL,'--','Color',[0.75 0.75 0.75],'LineWidth',1.2);
yline(ax2, base.IPL,'--','Color',[0.75 0.75 0.75],'LineWidth',1.2);
yline(ax3, base.GCL,'--','Color',[0.75 0.75 0.75],'LineWidth',1.2);

% ---- Shared y-axis for top row ----
yl1 = ylim(ax1); yl2 = ylim(ax2); yl3 = ylim(ax3);
ymin_top = min([yl1(1), yl2(1), yl3(1)]);
ymax_top = max([yl1(2), yl2(2), yl3(2)]);
yl_top   = [ymin_top, ymax_top];
yt_top   = ymin_top:20:ymax_top;
ylim(ax1, yl_top); yticks(ax1, yt_top);
ylim(ax2, yl_top); yticks(ax2, yt_top);
ylim(ax3, yl_top); yticks(ax3, yt_top);

%% ---------------------- Bottom row: split-dose comparison ----------------------
ax4 = axes('Parent',fig,'Position',[x1 botY axW axH]); hold(ax4,'on');
h_single_OPL = plot_layer_model(ax4, O_single, t_single, col.OPL, base.OPL);
h_6_OPL      = plot(ax4, t_6,  O_6,  '--', 'Color', col.OPL, 'LineWidth', 1.5);
h_12_OPL     = plot(ax4, t_12, O_12, ':',  'Color', col.OPL, 'LineWidth', 1.5);
h_24_OPL     = plot(ax4, t_24, O_24, '-.', 'Color', col.OPL, 'LineWidth', 1.5);
title(ax4,'Outer Plexiform Layer (OPL)','FontWeight','bold','FontSize',titleFS);

ax5 = axes('Parent',fig,'Position',[x2 botY axW axH]); hold(ax5,'on');
h_single_IPL = plot_layer_model(ax5, I_single, t_single, col.IPL, base.IPL);
h_6_IPL      = plot(ax5, t_6,  I_6,  '--', 'Color', col.IPL, 'LineWidth', 1.5);
h_12_IPL     = plot(ax5, t_12, I_12, ':',  'Color', col.IPL, 'LineWidth', 1.5);
h_24_IPL     = plot(ax5, t_24, I_24, '-.', 'Color', col.IPL, 'LineWidth', 1.5);
title(ax5,'Inner Plexiform Layer (IPL)','FontWeight','bold','FontSize',titleFS);

ax6 = axes('Parent',fig,'Position',[x3 botY axW axH]); hold(ax6,'on');
h_single_GCL = plot_layer_model(ax6, G_single, t_single, col.GCL, base.GCL);
h_6_GCL      = plot(ax6, t_6,  G_6,  '--', 'Color', col.GCL, 'LineWidth', 1.5);
h_12_GCL     = plot(ax6, t_12, G_12, ':',  'Color', col.GCL, 'LineWidth', 1.5);
h_24_GCL     = plot(ax6, t_24, G_24, '-.', 'Color', col.GCL, 'LineWidth', 1.5);
title(ax6,'Ganglion Cell Layer & Nerve Fiber Layer (GCL & NFL)','FontWeight','bold','FontSize',titleFS);

% ---- Shared y-axis for bottom row (with padding for legends) ----
yl4 = ylim(ax4); yl5 = ylim(ax5); yl6 = ylim(ax6);
ymin_bot = min([yl4(1), yl5(1), yl6(1)]);
ymax_bot = max([yl4(2), yl5(2), yl6(2)]);
yl_bot   = [ymin_bot, ymax_bot * 1.16];
yt_bot   = ymin_bot:20:ymax_bot;
ylim(ax4, yl_bot); yticks(ax4, yt_bot);
ylim(ax5, yl_bot); yticks(ax5, yt_bot);
ylim(ax6, yl_bot); yticks(ax6, yt_bot);

%% ---------------------- Legends (bottom row) ----------------------
legLabels = {'Single bolus (0 h)', 'Split dose (0 + 6 h)', ...
             'Split dose (0 + 12 h)', 'Split dose (0 + 24 h)'};

legBot4 = legend(ax4, [h_single_OPL, h_6_OPL, h_12_OPL, h_24_OPL], ...
    legLabels, 'Location','northeast','Box','off','FontSize',8);
legBot4.AutoUpdate = 'off';

legBot5 = legend(ax5, [h_single_IPL, h_6_IPL, h_12_IPL, h_24_IPL], ...
    legLabels, 'Location','northeast','Box','off','FontSize',8);
legBot5.AutoUpdate = 'off';

legBot6 = legend(ax6, [h_single_GCL, h_6_GCL, h_12_GCL, h_24_GCL], ...
    legLabels, 'Location','northeast','Box','off','FontSize',8);
legBot6.AutoUpdate = 'off';

%% ---------------------- Horizontal legend between rows ----------------------
axLegTop = axes('Parent',fig,'Position',[legX legY legW legH],'Visible','off');
hold(axLegTop,'on');
for k = 1:nMult
    hdlLegTop(k) = plot(axLegTop, NaN, NaN, 'Color', cmap(k,:), 'LineWidth', 1.8);
end
legTop = legend(axLegTop, hdlLegTop, lgdEntriesTop, ...
    'Orientation','horizontal', ...
    'NumColumns', nMult, ...
    'Location','north', ...
    'Box','off', ...
    'FontSize', 9.5, ...
    'FontWeight','bold');
legTop.AutoUpdate    = 'off';
legTop.ItemTokenSize = [13 9];

drawnow;

%% ====================== Local functions ======================
function [t_out, y_out] = simulate_split_peritoneum(p, IC_single, R0_total, tspan, t_inj2, opts_ode)
    IC1    = IC_single;
    IC1(4) = R0_total / 2;
    [t1, y1] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p), [tspan(1) t_inj2], IC1, opts_ode);
    y_inj    = y1(end,:);
    y_inj(4) = y_inj(4) + R0_total / 2;
    [t2, y2] = ode15s(@(tt,yy) rhs_microglia(tt,yy,p), [t_inj2 tspan(2)], y_inj, opts_ode);
    t_out = [t1; t2(2:end)];
    y_out = [y1; y2(2:end,:)];
end

function dydt = rhs_microglia(~, y, p)
    O = y(1); I = y(2); G = y(3); R = y(4); B = y(5);
    hill_G = B / (B + p.eta_G);
    hill_O = B / (B + p.eta_O);
    flux_G = p.k_G * G * (1 - G / p.G0);
    dO = -p.k_OI * O * hill_G + p.k_IO * I * hill_O;
    dI =  p.k_OI * O * hill_G - p.k_IO * I * hill_O - p.k_IG * I * hill_G - flux_G;
    dG =  p.k_IG * I * hill_G + flux_G;
    dR = -p.k_R * R;
    dB =  p.F * p.k_R * (p.Vol_P / p.Vol_B) * R - p.k_B * B;
    dydt = [dO; dI; dG; dR; dB];
end

function h_main = plot_layer_model(ax, y_sim, t_sim, colorRGB, baselineVal)
    h_main = plot(ax, t_sim, y_sim, 'Color', colorRGB, 'LineWidth', 2);
    if isfinite(baselineVal)
        yline(ax, baselineVal,'--','Color',[0.75 0.75 0.75],'LineWidth',1.2);
    end
    fmt_axis(ax);
    xlabel(ax,'Time (h)');
    ylabel(ax,'Microglia density (cells/mm^2)');
end

function fmt_axis(ax)
    set(ax,'LineWidth',2,'FontWeight','bold','TickDir','out');
    box(ax,'off');
    try, ax.Toolbar.Visible = 'off'; catch, end
end