%% Retina_IV_vs_IP_StylizedOverlay_InsetBlood.m
% Solid  = IV injection (biphasic blood model)
% Dashed = IP injection (peritoneum -> blood model)
%
% Layout (2x2):
%   Top-left     = Blood (IP main + IV inset)
%   Top-right    = OPL
%   Bottom-left  = IPL
%   Bottom-right = GCL/NFL
%
% CHANGE REQUEST:
%   Remove "Data" from legend in the inset panel.

clear; clc; close all;

%% ====================== Load data ======================
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

tdata = S.time(:); % hours

OPL_mean = mean(S.OPL_individual,2,'omitnan');
IPL_mean = mean(S.IPL_individual,2,'omitnan');
GCL_mean = mean(S.GCL_NFL_individual,2,'omitnan');

O0 = OPL_mean(1);
I0 = IPL_mean(1);
G0 = GCL_mean(1);

%% ====================== Shared microglia parameters ======================
p = struct();
p.k_OI  = 0.221192;
p.k_IO  = 0.42835;
p.k_IG  = 0.0572969;
p.k_G   = 0.0610402;

p.eta_O = 1667;
p.eta_G = 210;
p.G0    = G0;

%% ====================== IV BLOOD ======================
k1_h = 60 * 0.30;
k2_h = 60 * 0.03;

Btot0 = 150000;
B1_0  = 0.7 * Btot0;
B2_0  = 0.3 * Btot0;

paramsIV = struct();
paramsIV.micro = p;
paramsIV.k1_h  = k1_h;
paramsIV.k2_h  = k2_h;

IC_IV = [O0; I0; G0; B1_0; B2_0];

%% ====================== IP PARAMETERS ======================
paramsIP = p;
paramsIP.k_R   = 1.71;
paramsIP.F     = 0.03;
paramsIP.Vol_P = 0.06;
paramsIP.Vol_B = 0.9548;
paramsIP.k_B   = 0.05776;

IC_IP = [O0; I0; G0; 2021000; 0];

%% ====================== Solve models ======================
opts  = odeset('RelTol',1e-9,'AbsTol',1e-11);
t_end = max(tdata(end), 50);
tspan = [tdata(1) t_end];

solIV = ode15s(@(t,y) rhs_IV(t,y,paramsIV), tspan, IC_IV, opts);
solIP = ode15s(@(t,y) rhs_IP(t,y,paramsIP), tspan, IC_IP, opts);

tplot = linspace(tspan(1), tspan(2), 800);

Yiv = deval(solIV, tplot).';
Yip = deval(solIP, tplot).';

Oiv = Yiv(:,1); Iiv = Yiv(:,2); Giv = Yiv(:,3);
Oip = Yip(:,1); Iip = Yip(:,2); Gip = Yip(:,3);

tInset_h  = linspace(0, 50, 500);
YBip_full = deval(solIP, tInset_h).';
Bip_full  = YBip_full(:,5);

%% ====================== Inset fit ======================
data_fit = [
    0.40864158163264896, 99.60937499999999
    2.43558673469386800, 66.6015625
    4.94052933673469100, 48.046875
   14.99872448979591900, 26.757812499999996
   30.06744260204081700, 17.578124999999932
];

t_fit_min = data_fit(:,1);
y_fit_pct = data_fit(:,2);

model_fit = @(pp,t) pp(1).*exp(-pp(3).*t) + pp(2).*exp(-pp(4).*t);
p0_fit = [0.7*y_fit_pct(1), 0.3*y_fit_pct(1), 0.30, 0.03];

useLSQ = exist('lsqcurvefit','file') == 2;
if useLSQ
    optLSQ = optimoptions('lsqcurvefit','Display','off');
    pHat = lsqcurvefit(model_fit, p0_fit, t_fit_min, y_fit_pct);
else
    sse_q = @(q) sum((model_fit(exp(q), t_fit_min) - y_fit_pct).^2);
    qHat = fminsearch(sse_q, log(p0_fit));
    pHat = exp(qHat);
end

Ahat=pHat(1); Bhat=pHat(2); k1fit=pHat(3); k2fit=pHat(4);
if k2fit>k1fit
    [Ahat,Bhat]=deal(Bhat,Ahat);
    [k1fit,k2fit]=deal(k2fit,k1fit);
end

tt_fit = linspace(0, max(t_fit_min), 600).';
B1_fit = Ahat .* exp(-k1fit .* tt_fit);
B2_fit = Bhat .* exp(-k2fit .* tt_fit);
Btot_fit = B1_fit + B2_fit;

%% ====================== Plotting ======================
col.OPL=[52 166 213]/255;
col.IPL=[83 206 144]/255;
col.GCL=[94 71 181]/255;
col.B=[241 99 85]/255;

lw=2.4;

fig = figure('Color','w','Position',[120 120 900 650]);
tlo = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

%% ------------------ Blood ------------------
axBlood = nexttile(tlo,1); hold(axBlood,'on');
plot(axBlood, tInset_h, Bip_full, '--','Color',col.B,'LineWidth',lw);

xlim(axBlood,[0 50]);
title(axBlood,'Blood (B)','FontWeight','bold');
xlabel(axBlood,'Time (h)');
ylabel(axBlood,'LPS concentration (ng/mL)');
fmt(axBlood);

legend(axBlood,{'IP Injection'},'Location','northeast','Box','on');

%% ---- Inset ----
inset = axes('Position', inset_position(axBlood,[0.597 0.36 0.38 0.48]));
hold(inset,'on');

col_fast = col.B + (1-col.B)*0.45;
col_slow = col.B*0.75;

hTot  = plot(inset, tt_fit, Btot_fit,'-','Color',col.B,'LineWidth',2.4);

% Plot data but exclude from legend
plot(inset, t_fit_min, y_fit_pct,'o','Color',col.B,...
    'MarkerFaceColor','none','LineWidth',1.6);

hFast = plot(inset, tt_fit, B1_fit,'--','Color',col_fast,'LineWidth',2.0);
hSlow = plot(inset, tt_fit, B2_fit,':','Color',col_slow,'LineWidth',2.0);

xlim(inset,[0 30]);
xlabel(inset,'Time (min)','FontWeight','bold');
ylabel(inset,'Percent of Initial LPS','FontWeight','bold');
set(inset,'FontWeight','bold','LineWidth',1.5,'TickDir','out','Box','on');

legend(inset,[hTot hFast hSlow],...
    {'IV Injection','Fast','Slow'},...
    'Location','northeast','Box','on','FontWeight','bold','FontSize',8);

%% ------------------ OPL ------------------
axOPL = nexttile(tlo,2); hold(axOPL,'on');
plot(axOPL,tplot,Oiv,'-','Color',col.OPL,'LineWidth',lw);
plot(axOPL,tplot,Oip,'--','Color',col.OPL,'LineWidth',lw);
yline(axOPL,O0,'--','Color',[.75 .75 .75]);

title(axOPL,'Outer Plexiform Layer (OPL)','FontWeight','bold');
xlabel(axOPL,'Time (h)');
ylabel(axOPL,'Microglia density (cells/mm^2)');
fmt(axOPL);

legOPL = legend(axOPL,{'IV Injection','IP Injection','Baseline'}, ...
    'Location','east', 'Box','on');
set(legOPL,'FontWeight','bold');

%% ------------------ IPL ------------------
axIPL = nexttile(tlo,3); hold(axIPL,'on');
plot(axIPL,tplot,Iiv,'-','Color',col.IPL,'LineWidth',lw);
plot(axIPL,tplot,Iip,'--','Color',col.IPL,'LineWidth',lw);
yline(axIPL,I0,'--','Color',[.75 .75 .75]);

title(axIPL,'Inner Plexiform Layer (IPL)','FontWeight','bold');
xlabel(axIPL,'Time (h)');
ylabel(axIPL,'Microglia density (cells/mm^2)');
fmt(axIPL);

legIPL = legend(axIPL,{'IV Injection','IP Injection','Baseline'}, ...
    'Location','southeast', 'Box','on');
set(legIPL,'FontWeight','bold');

%% ------------------ GCL/NFL ------------------
axGCL = nexttile(tlo,4); hold(axGCL,'on');
plot(axGCL,tplot,Giv,'-','Color',col.GCL,'LineWidth',lw);
plot(axGCL,tplot,Gip,'--','Color',col.GCL,'LineWidth',lw);
yline(axGCL,G0,'--','Color',[.75 .75 .75]);

title(axGCL,'Ganglion Cell Layer / Nerve Fibre Layer (GCL/NFL)','FontWeight','bold');
xlabel(axGCL,'Time (h)');
ylabel(axGCL,'Microglia density (cells/mm^2)');
fmt(axGCL);

legGCL = legend(axGCL,{'IV Injection','IP Injection','Baseline'}, ...
    'Location','east', 'Box','on');
set(legGCL,'FontWeight','bold');

%% ====================== RHS ======================
function dydt = rhs_IV(~,y,P)
O=y(1); I=y(2); G=y(3); B1=y(4); B2=y(5);
B=B1+B2;
dB1=-P.k1_h*B1;
dB2=-P.k2_h*B2;
[dO,dI,dG]=micro_rhs(O,I,G,B,P.micro);
dydt=[dO;dI;dG;dB1;dB2];
end

function dydt = rhs_IP(~,y,p)
O=y(1); I=y(2); G=y(3); R=y(4); B=y(5);
[dO,dI,dG]=micro_rhs(O,I,G,B,p);
dR=-p.k_R*R;
dB=p.F*p.k_R*(p.Vol_P/p.Vol_B)*R - p.k_B*B;
dydt=[dO;dI;dG;dR;dB];
end

function [dO,dI,dG]=micro_rhs(O,I,G,B,p)
hillG=B/(B+p.eta_G);
hillO=B/(B+p.eta_O);
fluxG=p.k_G*G*(1-G/p.G0);

dO=-p.k_OI*O*hillG + p.k_IO*I*hillO;
dI=p.k_OI*O*hillG - p.k_IO*I*hillO - p.k_IG*I*hillG - fluxG;
dG=p.k_IG*I*hillG + fluxG;
end

function fmt(ax)
set(ax,'LineWidth',2,'FontWeight','bold','TickDir','out');
box(ax,'off');
end

function pos=inset_position(parentAx,rel)
pa=parentAx.Position;
pos=[pa(1)+rel(1)*pa(3), pa(2)+rel(2)*pa(4), ...
     rel(3)*pa(3), rel(4)*pa(4)];
end