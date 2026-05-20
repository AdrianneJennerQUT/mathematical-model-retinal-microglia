%% Paired PRCC-style bar plot
% Customizable bar heights for two conditions

clear; clc; close all;

%% ---------------------- USER INPUTS ------------------------------------
% Parameter labels (LaTeX allowed)
param_labels = { ...
    'k_{OI}', 'k_{IO}', 'k_{IG}', 'k_G', ...
    '\eta_O', '\eta_G', 'k_R', 'k_B', 'F', 'V_P', 'V_B'};

% Bar heights for set 1 (e.g. IPL AUC)
vals1 = [ ...
     0.7044,  -0.652497,  -0.20916,  0.327356,  0.590454, -0.389664, ...
     0.00566653,  0.620522,  -0.202995,  -0.235567,  0.212053];

% Bar heights for set 2 (e.g. IPL AUC, or GCL AUC)
vals2 = [ ...
     0.77716,  -0.741526,  -0.344499, 0.39878,  0.562454, -0.257293, ...
     0.0059392,  0.503991,  -0.300179, -0.352824,  0.310323];

% Titles / legend text
plot_title   = 'Sensitivity of IPL to Extended Model Parameters';
legend_names = {'IPL Density at 48 h','IPL AUC'};

%% ---------------------- BUILD GROUPED BAR DATA -------------------------
vals1 = vals1(:);
vals2 = vals2(:);

if numel(vals1) ~= numel(param_labels) || numel(vals2) ~= numel(param_labels)
    error('vals1, vals2, and param_labels must have the same length.');
end

Y = [vals1 vals2];    % N x 2 matrix for grouped bars

%% ---------------------- PLOT -------------------------------------------
fig = figure('Color','w','Position',[100 100 900 360]);

% Grouped bar chart
h = bar(Y, 'grouped', 'LineWidth', 1.5);
hold on;

col1 = [ 39 173  80]/255;  % darker green
col2 = [122 204 135]/255;  % lighter green

h(1).FaceColor = col1;
h(2).FaceColor = col2;

% Axis style
set(gca, 'XTick', 1:numel(param_labels), ...
         'XTickLabel', param_labels, ...
         'XTickLabelRotation', 45, ...
         'FontWeight', 'bold', ...
         'LineWidth', 1.5, ...
         'TickDir', 'out', ...
         'Box', 'off');

ylim([-1 1]);
yline(0, 'k-', 'LineWidth', 1.2);

ylabel('PRCC', 'FontWeight', 'bold');
title(plot_title, 'FontWeight', 'bold');

%% ---------------------- Legend inside plot, under title -------------------------
% Create legend first
leg = legend(h, legend_names, ...
    'Orientation','horizontal', ...
    'FontWeight','bold', ...
    'Box','off');

% Move legend into top-center of axes (inside plot area)
leg.Units = 'normalized';
leg.Position = [0.67, 0.87, 0.30, 0.05];  
% [x_center, y, width, height]

% Improve title spacing so legend fits nicely below it
title(plot_title, 'FontWeight','bold', 'Units','normalized', 'Position',[0.5, 1.07, 0]);

% Tighten axes to avoid shifting downward
set(gca, 'Position', [0.07 0.18 0.90 0.70]); 