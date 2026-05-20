%% Combined Figure: Row1 = IV Fit + Grouped Violin | Rows 2-4 = Region/KDE/Scenario
clear; clc; close all;
rng(1);

%% ====================== USER KNOBS ======================
nSamples  = 1000;
nKDE      = 400;
supp      = 'positive';
padSD     = 4;
nEval     = 800;
NBINS_END = 18;
LW_AX     = 2.0;
LW_BAR    = 1.26;
LW_KDE    = 2.1;
LW_SCEN   = 2.1;
LW_CURV   = 2.1;

%% ====================== LOAD DATA ======================
matfile = 'Retinal_data.mat';
if ~exist(matfile,'file')
    error('File "%s" not found in %s', matfile, pwd);
end
S        = load(matfile);
tdata    = S.time(:);
OPL_data = S.OPL_individual;
IPL_data = S.IPL_individual;
GCL_data = S.GCL_NFL_individual;

rowmean = @(A) mean(A,2,'omitnan');
baseOPL = rowmean(OPL_data); baseOPL = baseOPL(1);
baseIPL = rowmean(IPL_data); baseIPL = baseIPL(1);
baseGCL = rowmean(GCL_data); baseGCL = baseGCL(1);

tspan = [tdata(1) tdata(end)];
tEval = linspace(tspan(1), tspan(2), nEval);

%% ====================== IV BLOOD DATA ======================
iv_data = [
0.842696629213485,  82.46789727126806
1.9662921348314626, 66.53932584269663
3.089887640449433,  54.41860465116279
3.932584269662925,  48.83930464278143
4.494382022471903,  45.440958236167056
19.9438202247191,   26.40625
180.25503122831367, 19.71131158917417
227.34559333795977, 15.553435114503802
299.50345650419337, 16.6149068323
];
t_min = iv_data(:,1);
y_pct = iv_data(:,2);

iv_model = @(p,t) p(1).*exp(-p(3).*t) + p(2).*exp(-p(4).*t);
y0_iv    = y_pct(1);
p0_iv    = [0.7*y0_iv, 0.3*y0_iv, 0.30, 0.03];

useLSQ = exist('lsqcurvefit','file') == 2;
if useLSQ
    opts_lsq = optimoptions('lsqcurvefit','Display','off', ...
        'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);
    pHat = lsqcurvefit(iv_model, p0_iv, t_min, y_pct, ...
        [0,0,0,0],[Inf,Inf,Inf,Inf], opts_lsq);
else
    sse_q = @(q) sum((iv_model(exp(q),t_min)-y_pct).^2);
    qHat  = fminsearch(sse_q, log(max(p0_iv,1e-12)), optimset('Display','off'));
    pHat  = exp(qHat);
end
A_iv=pHat(1); B_iv=pHat(2); k1_iv=pHat(3); k2_iv=pHat(4);
if k2_iv > k1_iv
    [A_iv,B_iv]   = deal(B_iv,A_iv);
    [k1_iv,k2_iv] = deal(k2_iv,k1_iv);
end
yHat_iv  = iv_model([A_iv,B_iv,k1_iv,k2_iv], t_min);
res_iv   = y_pct - yHat_iv;
R2_iv    = 1 - sum(res_iv.^2)/sum((y_pct-mean(y_pct)).^2);
RMSE_iv  = sqrt(mean(res_iv.^2));
t12_fast = log(2)/k1_iv;
t12_slow = log(2)/k2_iv;
statsStr_short = sprintf('RMSE=%.3g\nk_1=%.4g /min  (t_{1/2}=%.1f)\nk_2=%.4f /min  (t_{1/2}=%.1f)', ...
    RMSE_iv, k1_iv, t12_fast, k2_iv, t12_slow);

tt_iv   = linspace(0,max(t_min),600).';
B1_iv   = A_iv.*exp(-k1_iv.*tt_iv);
B2_iv   = B_iv.*exp(-k2_iv.*tt_iv);
Btot_iv = B1_iv + B2_iv;

%% ====================== VIOLIN DATA ======================
GCLFovea  = [238.674;119.337;116.022;79.558;61.878;54.144;36.464;28.729;12.155;4.420;4.420;6.630];
IPLFovea  = [430.939;417.680;406.630;387.845;344.751;333.702;330.387;319.337;319.337;278.453;286.188];
OPLFovea  = [243.094;219.890;228.729;230.939;220.994;220.994;148.066;137.017;132.597;144.751;206.630;206.630];
GCLMacula = [197.790;156.906;103.867;98.343;103.867;65.193;54.144;57.459;43.094;32.044;39.779;46.409];
IPLMacula = [329.282;312.707;302.762;287.293;265.193;235.359;213.260;198.895;148.066;125.967;125.967;114.917];
OPLMacula = [155.801;151.381;162.431;134.807;133.702;133.702;124.862;103.867;93.923;82.873;98.343;86.188];
GCLPeri   = [101.657;108.287;98.343;90.608;96.133;92.818;79.558;57.459;79.558;76.243;71.823;79.558];
IPLPeri   = [282.873;250.829;239.779;220.994;195.580;195.580;173.481;137.017;141.436;137.017;125.967;122.652];
OPLPeri   = [191.160;181.215;170.166;166.851;170.166;173.481;148.066;151.381;93.923;86.188;90.608;97.238];

%% ====================== COLORS ======================
purple      = [241  99  85]/255;
purple_fast = purple + (1-purple)*0.45;
purple_slow = purple * 0.75;

colPeriphery      = [203 253 201]/255;
colFovea          = [  1 130   1]/255;
colMacula         = [  2 254   1]/255;
colMouseIV        = [0 0 0];
regionCurveColors = [colPeriphery; colFovea; colMacula];

colOPL      = [ 52 166 213]/255;
colIPL      = [ 83 206 144]/255;
colGCL      = [ 94  71 181]/255;
layerColors_kde = {colOPL; colIPL; colGCL};

plotColors = [59 50 92; 63 111 151; 85 175 163; 252 240 179]/255;
greyBL     = [0.75 0.75 0.75];

regionNames = {'Fovea','Macula','Periphery'};
layerNames  = {'OPL','IPL','GCL'};
layerTitles = {'OPL','IPL','GCL'};

%% ====================== GROUPED VIOLIN SETUP ======================
groupData = {
    {GCLFovea,  IPLFovea,  OPLFovea};
    {GCLMacula, IPLMacula, OPLMacula};
    {GCLPeri,   IPLPeri,   OPLPeri}
};
groupNames_v  = {'Fovea','Macula','Periphery'};
layerNames_v  = {'GCL/NFL','IPL','OPL'};
layerColors_v = {colGCL; colIPL; colOPL};

vSpacing = 0.8;
vGroupSpacing = 0.4;
nGroups = 3;
nLayers = 3;
vCentres = (0:nGroups-1) * (nLayers*vSpacing + vGroupSpacing);
vOffsets = (-(nLayers-1)/2 : 1 : (nLayers-1)/2) * vSpacing;
vxPos = zeros(nGroups, nLayers);
for g = 1:nGroups
    vxPos(g,:) = vCentres(g) + vOffsets;
end

%% ====================== SHARED MICROGLIA PARAMETERS ======================
p_shared       = struct();
p_shared.eta_O = 1667;   p_shared.eta_G = 210;
p_shared.k_OI  = 0.221192; p_shared.k_IO = 0.42835;
p_shared.k_IG  = 0.0572969; p_shared.k_G = 0.0610402;

A_amp=71.45873285588064; B_amp=27.106559432654475;
k1_min=0.3027261887787264; k2_min=0.0019131529006041971;
k1_h_fit=k1_min*60; k2_h_fit=k2_min*60;
B0_monkey=174600;
frac1=A_amp/(A_amp+B_amp); frac2=B_amp/(A_amp+B_amp);
B1_0_fit=frac1*B0_monkey; B2_0_fit=frac2*B0_monkey;

%% ====================== ROW 2 IC PARAMETERS ======================
IC_regions(1).name='Periphery'; IC_regions(1).O0=142;   IC_regions(1).I0=184;    IC_regions(1).G0=83;
IC_regions(2).name='Fovea';     IC_regions(2).O0=193.0; IC_regions(2).I0=349.9;  IC_regions(2).G0=62.3;
IC_regions(3).name='Macula';    IC_regions(3).O0=117.6; IC_regions(3).I0=221.08; IC_regions(3).G0=81.3;
nR=numel(IC_regions);

O0_mouse=mean(OPLPeri,'omitnan'); I0_mouse=mean(IPLPeri,'omitnan'); G0_mouse=mean(GCLPeri,'omitnan');
k1_h_PK=60*0.30; k2_h_PK=60*0.03;
Btot0_PK=2.4889e5*9; B1_0_PK=0.7*Btot0_PK; B2_0_PK=0.3*Btot0_PK;
pPK=p_shared; pPK.G0=G0_mouse; pPK.k1_h=k1_h_PK; pPK.k2_h=k2_h_PK;
IC_PK=[O0_mouse;I0_mouse;G0_mouse;B1_0_PK;B2_0_PK];

paramsFit=struct(); paramsFit.micro=p_shared;
paramsFit.micro.k1_h=k1_h_fit; paramsFit.micro.k2_h=k2_h_fit;
paramsFit.k1_h=k1_h_fit; paramsFit.k2_h=k2_h_fit;

%% ====================== ROW 4 SCENARIO PARAMETERS ======================
p_sc=struct(); p_sc.k_R=1.71; p_sc.F=0.03;
p_sc.Vol_P=0.06; p_sc.Vol_B=0.9548; p_sc.k_B=0.05776;
p_sc.eta_O=1667; p_sc.eta_G=210; p_sc.k_OI=0.221192;
p_sc.k_IO=0.42835; p_sc.k_IG=0.0572969; p_sc.k_G=0.0610402;

O0_base=baseOPL; I0_base=baseIPL; G0_base=baseGCL;
Mtot=O0_base+I0_base+G0_base;
R0=2.4889e+05*9; B0_sc=0;
fO_base=O0_base/Mtot; fI_base=I0_base/Mtot;
dom=0.66; rem_frac=(1-dom)/2;
sc_names={'Homeostatic','OPL dominant','IPL dominant','GCL dominant'};
sc_fO=[fO_base,dom,rem_frac,rem_frac];
sc_fI=[fI_base,rem_frac,dom,rem_frac];
nSc=numel(sc_names);

%% ====================== COMPUTATION: ROW 2 ======================
opts_ode=odeset('RelTol',1e-9,'AbsTol',1e-11,'OutputFcn',[],'Stats','off', ...
    'MaxStep',(tspan(2)-tspan(1))/500);
tgrid=linspace(tspan(1),tspan(2),800);
t_all=cell(nR,1); O_all=cell(nR,1); I_all=cell(nR,1); G_all=cell(nR,1);
for k=1:nR
    pk=p_shared; pk.G0=IC_regions(k).G0; pk.k1_h=k1_h_fit; pk.k2_h=k2_h_fit;
    IC_reg=[IC_regions(k).O0;IC_regions(k).I0;IC_regions(k).G0;B1_0_fit;B2_0_fit];
    [t,y]=ode15s(@(tt,yy) rhs_microglia(tt,yy,pk),tgrid,IC_reg,opts_ode);
    t_all{k}=t; O_all{k}=y(:,1); I_all{k}=y(:,2); G_all{k}=y(:,3);
end
[t_PK,y_PK]=ode15s(@(tt,yy) rhs_microglia(tt,yy,pPK),tgrid,IC_PK,opts_ode);
O_PK=y_PK(:,1); I_PK=y_PK(:,2); G_PK=y_PK(:,3);
fprintf('Row 2 done.\n');

%% ====================== COMPUTATION: ROW 3 ======================
opts_mc=odeset('RelTol',1e-8,'AbsTol',1e-10,'OutputFcn',[],'Stats','off');
allData={OPLFovea,IPLFovea,GCLFovea; OPLMacula,IPLMacula,GCLMacula; OPLPeri,IPLPeri,GCLPeri};
O_end_all=cell(3,1); I_end_all=cell(3,1); G_end_all=cell(3,1);
for reg=1:3
    OPL_r=allData{reg,1}(:); IPL_r=allData{reg,2}(:); GCL_r=allData{reg,3}(:);
    OPL_samp=max(0,sample_from_kde(OPL_r,nSamples,nKDE,supp,padSD));
    IPL_samp=max(0,sample_from_kde(IPL_r,nSamples,nKDE,supp,padSD));
    GCL_samp=max(0,sample_from_kde(GCL_r,nSamples,nKDE,supp,padSD));
    O_end=nan(nSamples,1); I_end=nan(nSamples,1); G_end=nan(nSamples,1);
    for s=1:nSamples
        O0=OPL_samp(s); I0=IPL_samp(s); G0=GCL_samp(s);
        paramsFit.micro.G0=G0;
        IC=[O0;I0;G0;B1_0_fit;B2_0_fit];
        sol=ode15s(@(t,y) rhs_IV(t,y,paramsFit),tspan,IC,opts_mc);
        Y=deval(sol,tEval).';
        O_end(s)=Y(end,1); I_end(s)=Y(end,2); G_end(s)=Y(end,3);
    end
    O_end_all{reg}=O_end; I_end_all{reg}=I_end; G_end_all{reg}=G_end;
    fprintf('KDE region %s done.\n',regionNames{reg});
end

%% ====================== COMPUTATION: ROW 4 ======================
opts_sc=odeset('RelTol',1e-8,'AbsTol',1e-10,'OutputFcn',[],'Stats','off');
tgrid_sc=linspace(tspan(1),tspan(2),800);
sc_t=cell(nSc,1); sc_O=cell(nSc,1); sc_I=cell(nSc,1); sc_G=cell(nSc,1);
for k=1:nSc
    fO=sc_fO(k); fI=sc_fI(k); fG=max(0,1-fO-fI);
    O0=fO*Mtot; I0=fI*Mtot; G0=fG*Mtot;
    IC=[O0 I0 G0 R0 B0_sc]; pk=p_sc; pk.G0=G0;
    [t,y]=ode15s(@(tt,yy) rhs_micro(tt,yy,pk),tgrid_sc,IC,opts_sc);
    sc_t{k}=t; sc_O{k}=y(:,1); sc_I{k}=y(:,2); sc_G{k}=y(:,3);
end
fprintf('All computation done. Building figure...\n');

%% ====================== FIGURE LAYOUT ======================
scale     =  0.78;   % reduce below 1.0 to shrink everything
leftMarg  =  round(85  * scale);
rightMarg =  round(30  * scale);
topMarg   =  round(50  * scale);
botMarg   =  round(75  * scale);
gapX      =  round(45  * scale);
gapR1R2   =  round(60  * scale);
gapRows   =  round(90  * scale);

nCols234     = 3;
panelW       = 260*0.65;   % increased — makes rows 2-4 wider
row234TotalW = nCols234*panelW + (nCols234-1)*gapX;
panelH       = panelW * (2/3);

% Row 1 sized independently
row1TotalW   = (260*0.65)*3 + 2*gapX;   % fixed at original width
panelH_1     = row1TotalW * 0.35;        % violin panel height
panelH_1_iv  = row1TotalW * 0.325;        % IV panel height — adjust independently
row1PanelW   = row1TotalW * 0.35;
gapR1        = 50;
violinPanelW = row1TotalW - row1PanelW - gapR1;
violinPanelX = leftMarg + row1PanelW + gapR1;

figW = leftMarg + row234TotalW + rightMarg;
figH = topMarg + panelH_1 + gapR1R2 + 3*panelH + 2*gapRows + botMarg;

row4B = botMarg;
row3B = row4B + panelH  + gapRows;
row2B = row3B + panelH  + gapRows;
row1B = row2B + panelH  + gapR1R2;

ivOffset = 14;   % pixels to shift IV panel upward — increase to move higher
pxR1_iv = [leftMarg, row1B + ivOffset, row1PanelW, panelH_1_iv];
pxR1_vio = [violinPanelX, row1B, violinPanelW, panelH_1];
pxR  = @(row,col) [leftMarg + (col-1)*(panelW+gapX), row, panelW, panelH];

set(0,'DefaultFigureVisible','off');
fig = figure('Color','w','Position',[60 60 figW figH]);

%% ====================== ROW 1, PANEL 1: IV FIT ======================
ax_iv = axes(fig,'Units','pixels','Position',pxR1_iv);
hold(ax_iv,'on');

hData_iv = plot(ax_iv,t_min,y_pct,'o','Color',purple, ...
    'MarkerEdgeColor',purple,'MarkerFaceColor','none','LineWidth',1.8*0.6,'MarkerSize',7*0.6);
hTot_iv  = plot(ax_iv,tt_iv,Btot_iv,'-' ,'Color',purple,     'LineWidth',3.0*0.7);

xlabel(ax_iv,'Time (min)','FontWeight','bold','FontSize',7.2);
set(ax_iv,'LineWidth',2,'FontWeight','bold','TickDir','out','FontSize',8);
box(ax_iv,'off');
text(ax_iv, ax_iv.XLim(1)-30, mean(ax_iv.YLim), 'Normalized LPS Concentration', ...
    'FontSize',7.2,'FontWeight','bold','Color','k', ...
    'HorizontalAlignment','center','Rotation',90,'VerticalAlignment','bottom');

universal=6.3;

legLabels_iv={'Data','Model'};
leg_iv=legend(ax_iv,[hData_iv,hTot_iv],legLabels_iv, ...
    'Location','northeast','Box','on','FontSize',universal,'FontWeight','bold');
leg_iv.ItemTokenSize = [7, 3.5];
leg_iv.LineWidth = 1;
drawnow;

leg_iv.Units = 'pixels';
legSz = leg_iv.Position(3:4);
leg_iv_X = pxR1_iv(1) + 13.78;                    % <-- move legend left/right
leg_iv_Y = row1B + panelH_1 - legSz(2) - 40;   % <-- move legend up/down
leg_iv.Position = [leg_iv_X, leg_iv_Y, legSz(1), legSz(2)];
drawnow;

stats_X = pxR1_iv(1) + 13;   % <-- move stats box left/right
stats_Y = row1B + 115;          % <-- move stats box up/down
stats_W = 120;                % <-- stats box width  (px)
stats_H = 44;                 % <-- stats box height (px)
statsBox_n = [stats_X, stats_Y, stats_W, stats_H] ./ [figW figH figW figH];
annotation(fig,'textbox',statsBox_n, ...
    'String',statsStr_short,'Interpreter','tex','FontWeight','bold', ...
    'FontSize',universal,'EdgeColor','k','LineWidth',1, ...
    'BackgroundColor','w','VerticalAlignment','top','Margin',2.5, ...
    'HorizontalAlignment','center');

monkey    = imread('thisisthemonkey.png');
axMonkey  = axes(fig,'Units','pixels','Position', ...
    [pxR1_iv(1)+pxR1_iv(3)*0.55*1.06, row1B+panelH_1*0.46*1.06, ...
     pxR1_iv(3)*0.46*0.90, panelH_1*0.68*0.90]);
imshow(monkey,'Parent',axMonkey);
axMonkey.Visible='off';
drawnow;

%% ====================== ROW 1, PANEL 2: GROUPED VIOLIN ======================
ax_v = axes(fig,'Units','pixels','Position',pxR1_vio);
hold(ax_v,'on'); box(ax_v,'off');
set(ax_v,'FontSize',9,'FontWeight','bold','LineWidth',1.8,'TickDir','out');
ax_v.YAxis.Visible = 'off';
ax_v.XAxis.Visible = 'off';

yMax_v  = 500;
yLims_v = [-25, 560];
ylim(ax_v, yLims_v);

vMargin = 0.35 + 0.5;
vxLeft  = min(vxPos(:)) - vMargin + 0.2;
vxRight = max(vxPos(:)) + vMargin - 0.2;
xlim(ax_v, [vxLeft, vxRight]);

% Grid lines
vGridGrey = [0.6 0.6 0.6];
vAxisCol  = [0.2 0.2 0.2];
yTicks_v  = 0:100:500;
for k = 1:numel(yTicks_v)
    plot(ax_v,[vxLeft vxRight],[yTicks_v(k) yTicks_v(k)],'-','Color',vGridGrey,'LineWidth',0.7);
end
uistack(findobj(ax_v,'Type','line'),'bottom');
plot(ax_v,[vxLeft vxRight],[0 0],'Color','k','LineWidth',1.4);

% Y spine + ticks
plot(ax_v,[vxLeft vxLeft],[0 yMax_v],'Color',vAxisCol,'LineWidth',1.4);
for k = 1:numel(yTicks_v)
    plot(ax_v,[vxLeft, vxLeft-0.12],[yTicks_v(k) yTicks_v(k)],'Color',vAxisCol,'LineWidth',1.2);
    text(ax_v, vxLeft-0.18, yTicks_v(k), sprintf('%d',yTicks_v(k)), ...
        'FontSize',7.2,'FontWeight','bold','Color','k', ...
        'HorizontalAlignment','right','VerticalAlignment','middle');
end
text(ax_v, vxLeft-0.64, mean([0 yMax_v]), 'Microglia density (cells/mm^{2})', ...
    'FontSize',7,'FontWeight','bold','Color','k', ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',90);
% Draw violins
vViolinWidth  = 0.35;
vViolinAlpha  = 0.72;
vViolinEdgeLW = 1.2;
vBoxWidth     = 0.10;
vBoxLW        = 1.8;
vMedianLW     = 2.8;
vWhiskLW      = 1.4;
vCapWidth     = 0.12;
vCapLW        = 1.6;
vPointSize    = 4;
vRainJitter   = 0.055;
vNKDE         = 512;
vKdeTaper     = 0.001;

hLeg_v = gobjects(nLayers,1);

for g = 1:nGroups
    for L = 1:nLayers
        x  = groupData{g}{L};
        x  = x(isfinite(x));
        x0 = vxPos(g,L);
        c  = layerColors_v{L};

        sd = std(x,'omitnan');
        if sd<=0||~isfinite(sd), sd=1; end
        lo = max(0, min(x)-3*sd); hi = max(x)+3*sd;
        yi = linspace(lo,hi,vNKDE);
        f  = ksdensity(x,yi,'Bandwidth',sd*0.4); f=f(:)';
        if max(f)>0, f=f./max(f); end
        keep=f>=vKdeTaper; yi=yi(keep); halfW=f(keep)*vViolinWidth;
        keep2=yi<=yMax_v; yi=yi(keep2); halfW=halfW(keep2);

        hv = fill(ax_v,[x0+halfW,fliplr(x0-halfW)],[yi,fliplr(yi)], ...
            c,'FaceAlpha',vViolinAlpha,'EdgeColor',c*0.6,'LineWidth',vViolinEdgeLW);
        if g==1, hLeg_v(L)=hv; end

        q=quantile(x,[0.25 0.5 0.75]); q1=q(1); med=q(2); q3=q(3);
        ylo=min(x); yhi=min(max(x),yMax_v);
        rectangle(ax_v,'Position',[x0-vBoxWidth/2,q1,vBoxWidth,q3-q1], ...
            'EdgeColor',c*0.55,'LineWidth',vBoxLW,'FaceColor',[1 1 1]*0.95);
        plot(ax_v,[x0-vBoxWidth/2,x0+vBoxWidth/2],[med med],'Color','k','LineWidth',vMedianLW);
        plot(ax_v,[x0 x0],[ylo q1],'Color',c*0.55,'LineWidth',vWhiskLW);
        plot(ax_v,[x0 x0],[q3 yhi],'Color',c*0.55,'LineWidth',vWhiskLW);
        plot(ax_v,[x0-vCapWidth/2,x0+vCapWidth/2],[ylo ylo],'Color',c*0.55,'LineWidth',vCapLW);
        plot(ax_v,[x0-vCapWidth/2,x0+vCapWidth/2],[yhi yhi],'Color',c*0.55,'LineWidth',vCapLW);
        jitter=(rand(size(x))-0.5)*2*vRainJitter;
        scatter(ax_v,x0+jitter,x,vPointSize, ...
            'MarkerFaceColor',c*0.65,'MarkerEdgeColor','w', ...
            'LineWidth',0.5,'MarkerFaceAlpha',0.85);
    end

    text(ax_v, mean(vxPos(g,:)), yLims_v(1)-8, groupNames_v{g}, ...
        'FontSize',8,'FontWeight','bold','Color','k', ...
        'HorizontalAlignment','center','VerticalAlignment','top');

    if g < nGroups
        xSep = mean([vxPos(g,end), vxPos(g+1,1)]);
        plot(ax_v,[xSep xSep],[0 yMax_v],'--','Color',[0.82 0.82 0.82],'LineWidth',0.8);
    end
end

leg_v = legend(ax_v, hLeg_v, layerNames_v, 'Box','off', ...
    'Orientation','horizontal','Location','northoutside');
set(leg_v,'FontSize',7,'FontWeight','bold');
drawnow;

% Manually pull legend closer to plot
leg_v.Units = 'pixels';
legSz_v = leg_v.Position(3:4);
leg_v.Position = [pxR1_vio(1) + (violinPanelW-legSz_v(1))/2, ...
                  row1B + panelH_1 - legSz_v(2) + 6, ...
                  legSz_v(1), legSz_v(2)];
drawnow;

%% ====================== ROW 2: REGION CURVES ======================
hLegRow2=gobjects(nR+1,1);
for col=1:3
    ax=axes(fig,'Units','pixels','Position',pxR(row2B,col)); %#ok<LAXES>
    if col==1, ax_row2_col1=ax; end
    hold(ax,'on');
    for k=1:nR
        switch col; case 1,ydata=O_all{k}; case 2,ydata=I_all{k}; case 3,ydata=G_all{k}; end
        h=plot(ax,t_all{k},ydata,'LineWidth',LW_CURV,'Color',regionCurveColors(k,:));
        if col==1, hLegRow2(k)=h; end
    end
    switch col; case 1,yPK=O_PK; case 2,yPK=I_PK; case 3,yPK=G_PK; end
    h=plot(ax,t_PK,yPK,'--','LineWidth',LW_CURV,'Color',colMouseIV);
    if col==1, hLegRow2(nR+1)=h; end
    title(ax,layerTitles{col},'FontWeight','bold');
    xlabel(ax,'Time (h)','FontWeight','bold','FontSize',7.2);
    if col==1, ylabel(ax,'Microglia density (cells/mm^2)','FontWeight','bold','FontSize',7.2); end
    xlim(ax,[0 48]); fmt_axis(ax,LW_AX);
end
legRow2Names={IC_regions.name,'Mouse IV'};
legRow2=legend(ax_row2_col1,hLegRow2,legRow2Names,'Orientation','horizontal','Box','off');
set(legRow2,'FontSize',7,'FontWeight','bold','AutoUpdate','off');
legRow2.Units='pixels';
legRow2.Position(1:2)=[leftMarg+105, row2B-45];

%% ====================== ROW 3: KDE HISTOGRAMS ======================
regionEndData={{O_end_all{1},I_end_all{1},G_end_all{1}}; ...
               {O_end_all{2},I_end_all{2},G_end_all{2}}; ...
               {O_end_all{3},I_end_all{3},G_end_all{3}}};
hKDE=gobjects(3,1);
for reg=1:3
    ax=axes(fig,'Units','pixels','Position',pxR(row3B,reg)); %#ok<LAXES>
    hold(ax,'on');
    for L=1:3
        x=regionEndData{reg}{L};
        lcol=layerColors_kde{L}; lcol_light=0.55*lcol+0.45*[1 1 1];
        hh=histogram(ax,x,NBINS_END,'Normalization','pdf','FaceColor',lcol_light, ...
            'EdgeColor',lcol,'LineWidth',LW_BAR,'FaceAlpha',0.55);
        hh.Annotation.LegendInformation.IconDisplayStyle='off';
        hKDE(L)=overlay_kde_handle(ax,x,nKDE,supp,padSD,lcol,LW_KDE);
    end
    title(ax,regionNames{reg},'FontWeight','bold');
    xlabel(ax,{'Cell Density at 48h (cells/mm^2)'},'FontWeight','bold','FontSize',7.2);
    if reg==1, ylabel(ax,'PDF','FontWeight','bold','FontSize',7.2); end
    lg=legend(ax,hKDE,layerNames,'Box','off','Location','northeast');
    set(lg,'FontSize',7.2,'FontWeight','bold');
    fmt_axis(ax,LW_AX); force_absolute_y(ax);
end

%% ====================== ROW 4: SCENARIO CURVES ======================
hLegRow4=gobjects(nSc,1);
scData={sc_O,sc_I,sc_G}; baseLines=[baseOPL,baseIPL,baseGCL];
for col=1:3
    ax=axes(fig,'Units','pixels','Position',pxR(row4B,col)); %#ok<LAXES>
    hold(ax,'on');
    for k=1:nSc
        h=plot(ax,sc_t{k},scData{col}{k},'LineWidth',LW_SCEN,'Color',plotColors(k,:));
        if col==1, hLegRow4(k)=h; end
    end
    yline(ax,baseLines(col),'--','Color',greyBL,'LineWidth',2);
    title(ax,layerTitles{col},'FontWeight','bold');
    xlabel(ax,'Time (h)','FontWeight','bold');
    if col==1, ylabel(ax,'Microglia density (cells/mm^2)','FontWeight','bold','FontSize',7.2); end
    xlim(ax,[0 48]); fmt_axis(ax,LW_AX);
end
legRow4=legend(hLegRow4,sc_names,'Orientation','horizontal','Box','off');
set(legRow4,'FontSize',7,'FontWeight','bold','AutoUpdate','off');
drawnow;
legRow4.Units='pixels';
legSz_r4 = legRow4.Position(3:4);
legRow4.Position(1:2) = [leftMarg + (row234TotalW - legSz_r4(1))/2, row4B-45];

%% ====================== LOCAL FUNCTIONS ======================
function xSamp=sample_from_kde(x,nSamp,nKDE,supp,padSD)
    x=x(isfinite(x)); mu=mean(x,'omitnan'); sd=std(x,'omitnan');
    if sd<=0||~isfinite(sd), sd=1; end
    xLo=min(max(0,mu-padSD*sd),min(x)); xHi=max(mu+padSD*sd,max(x));
    if ~(xHi>xLo), xSamp=repmat(mu,nSamp,1); return; end
    xi=linspace(xLo,xHi,nKDE);
    f=max(ksdensity(x,xi,'Support',supp),0);
    if sum(f)<=0, xSamp=x(randi(numel(x),nSamp,1)); return; end
    dx=xi(2)-xi(1); cdf=cumsum(f)*dx; cdf=cdf/cdf(end);
    [cdfU,ia]=unique(cdf,'stable'); xiU=xi(ia);
    cdfU(1)=0; cdfU(end)=1;
    xSamp=interp1(cdfU,xiU,rand(nSamp,1),'pchip');
    xSamp=max(min(xSamp,xiU(end)),xiU(1));
end

function h=overlay_kde_handle(ax,x,nKDE,supp,padSD,lineColor,lw)
    x=x(isfinite(x)); mu=mean(x,'omitnan'); sd=std(x,'omitnan');
    if sd<=0||~isfinite(sd), sd=1; end
    xLo=min(max(0,mu-padSD*sd),min(x)); xHi=max(mu+padSD*sd,max(x));
    xi=linspace(xLo,xHi,nKDE);
    f=ksdensity(x,xi,'Support',supp);
    h=plot(ax,xi,f,'-','Color',lineColor,'LineWidth',lw);
end

function dydt=rhs_microglia(~,y,p)
    O=y(1);I=y(2);G=y(3);B1=y(4);B2=y(5); B=B1+B2;
    hG=B./(B+p.eta_G); hO=B./(B+p.eta_O);
    fluxG=p.k_G*G.*(1-G./p.G0);
    dO=-p.k_OI*O.*hG+p.k_IO*I.*hO;
    dI= p.k_OI*O.*hG-p.k_IO*I.*hO-p.k_IG*I.*hG-fluxG;
    dG= p.k_IG*I.*hG+fluxG;
    dB1=-p.k1_h*B1; dB2=-p.k2_h*B2;
    dydt=[dO;dI;dG;dB1;dB2];
end

function dydt=rhs_IV(~,y,P)
    O=y(1);I=y(2);G=y(3);B1=y(4);B2=y(5); B=B1+B2;
    dB1=-P.k1_h*B1; dB2=-P.k2_h*B2;
    [dO,dI,dG]=micro_rhs(O,I,G,B,P.micro);
    dydt=[dO;dI;dG;dB1;dB2];
end

function [dO,dI,dG]=micro_rhs(O,I,G,B,p)
    hillG=B/(B+p.eta_G); hillO=B/(B+p.eta_O);
    fluxG=p.k_G*G*(1-G/p.G0);
    dO=-p.k_OI*O*hillG+p.k_IO*I*hillO;
    dI= p.k_OI*O*hillG-p.k_IO*I*hillO-p.k_IG*I*hillG-fluxG;
    dG= p.k_IG*I*hillG+fluxG;
end

function dydt=rhs_micro(~,y,p)
    O=y(1);I=y(2);G=y(3);R=y(4);B=y(5);
    hG=B/(B+p.eta_G); hO=B/(B+p.eta_O);
    flux=p.k_G*G*(1-G/p.G0);
    dO=-p.k_OI*O*hG+p.k_IO*I*hO;
    dI= p.k_OI*O*hG-p.k_IO*I*hO-p.k_IG*I*hG-flux;
    dG= p.k_IG*I*hG+flux;
    dR=-p.k_R*R; dB=p.F*p.k_R*(p.Vol_P/p.Vol_B)*R-p.k_B*B;
    dydt=[dO;dI;dG;dR;dB];
end

function fmt_axis(ax,LW)
    set(ax,'LineWidth',LW,'FontWeight','bold','TickDir','out','FontSize',7.2);
    box(ax,'off');
    try ax.Toolbar.Visible='off'; catch, end
end

function force_absolute_y(ax)
    ax.YAxis.Exponent=0;
end

%% ====================== EXPORT ======================
set(0,'DefaultFigureVisible','on');
fig.Visible='on'; figure(fig); drawnow;
print(fig,'combined_4row_figure','-dpng','-r150');
fprintf('Figure saved.\n');