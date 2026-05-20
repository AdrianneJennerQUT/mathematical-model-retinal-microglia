%% pk_fit_two_doses_log_space_with_stats.m
% Fits the PK model using log-transformed concentrations.
% Objective minimized:
%     log(model concentration) - log(observed concentration)

clear; clc;

%% ---------------------- Data ----------------------
time_points = [0.5, 2, 4, 12, 24]';

dose_0_5 = [136.986, 223.287, 206.849, 98.630, 60.274]';
dose_1_0 = [250.685, 372.603, 363.014, 242.466, 139.726]';

combined_dose_data = [dose_0_5; dose_1_0];

% Error bars for 0.5 mg/kg
high_points_05mg = [169.15, 260.64, 222.34, 119.15, 78.72]';
low_points_05mg  = [108.51, 192.55, 190.43, 81.91, 41.49]';

err05_lower = dose_0_5 - low_points_05mg;
err05_upper = high_points_05mg - dose_0_5;

% Error bars for 1.0 mg/kg
high_points_1mg = [274.47, 398.94, 394.68, 282.98, 152.13]';
low_points_1mg  = [219.15, 350, 345.74, 202.13, 125.53]';

err10_lower = dose_1_0 - low_points_1mg;
err10_upper = high_points_1mg - dose_1_0;

%% ---------------------- Fixed parameters ----------------------
Vol_P = 0.06;
Vol_B = 0.9548;
k2    = 0.05776;

%% ---------------------- Colors ----------------------
col_data_05  = [1, 0.5, 0.5];
col_data_10  = [0.5, 0.2, 0.9];

col_model_05 = [1, 0.75, 0.8];
col_model_10 = [0.7, 0.4, 0.9];

ci_alpha = 0.30;

%% ---------------------- Fit in log space ----------------------
params0 = [1, 0.05];  % [kPB, F]

% Log-space residuals:
% residual = log(model concentration) - log(observed concentration)
obj = @(p) log(FObjective_data_space(p, time_points, Vol_P, Vol_B, k2)) ...
         - log(combined_dose_data);

use_lsqnonlin = exist('lsqnonlin', 'file') == 2;

if use_lsqnonlin
    opts = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'FiniteDifferenceType', 'central', ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);

    [popt, ~, residual_log, ~, ~, ~, J] = lsqnonlin(obj, params0, [], [], opts);
else
    warning('lsqnonlin not found. Falling back to fminsearch.');
    opts = optimset('Display', 'iter', 'TolX', 1e-10, 'TolFun', 1e-10);

    popt = fminsearch(@(p) sum(obj(p).^2), params0, opts);

    residual_log = obj(popt);
    J = finiteDiffJacobian(obj, popt);
end

kPB_hat = popt(1);
F_hat   = popt(2);

%% ---------------------- Parameter confidence intervals ----------------------
N = numel(residual_log);
p = numel(popt);

SSE_log = sum(residual_log.^2);
s2 = SSE_log / max(N - p, 1);

pcov = s2 * inv(J.' * J);
perr = sqrt(diag(pcov));

z = 1.96;

kPB_lower = kPB_hat - z * perr(1);
kPB_upper = kPB_hat + z * perr(1);

F_lower = F_hat - z * perr(2);
F_upper = F_hat + z * perr(2);

fprintf('\nFitted parameters and 95%% CI from log-space fit:\n');
fprintf('k_PB = %.6g [%.6g, %.6g]\n', kPB_hat, kPB_lower, kPB_upper);
fprintf('F    = %.6g [%.6g, %.6g]\n', F_hat, F_lower, F_upper);

%% ---------------------- Goodness-of-fit statistics ----------------------
% RMSE and R^2 are evaluated on the original concentration scale.
model_0_5_obs = getModelOutput_matlab([kPB_hat, F_hat], time_points, 0.5, Vol_P, Vol_B, k2);
model_1_0_obs = getModelOutput_matlab([kPB_hat, F_hat], time_points, 1.0, Vol_P, Vol_B, k2);

model_all = [model_0_5_obs; model_1_0_obs];
data_all  = combined_dose_data;

residual_data = model_all - data_all;

RMSE_data = sqrt(mean(residual_data.^2));

SS_res = sum(residual_data.^2);
SS_tot = sum((data_all - mean(data_all)).^2);

R_squared = 1 - SS_res / SS_tot;

fprintf('\nGoodness-of-fit statistics evaluated on original concentration scale:\n');
fprintf('RMSE = %.6g ng/mL\n', RMSE_data);
fprintf('R^2  = %.6f\n', R_squared);

%% ---------------------- Plotting grid ----------------------
vector_time_points = linspace(0, 24, 1000)';

model_0_5 = getModelOutput_matlab([kPB_hat, F_hat], vector_time_points, 0.5, Vol_P, Vol_B, k2);
model_1_0 = getModelOutput_matlab([kPB_hat, F_hat], vector_time_points, 1.0, Vol_P, Vol_B, k2);

model_0_5_upper = getModelOutput_matlab([kPB_upper, F_upper], vector_time_points, 0.5, Vol_P, Vol_B, k2);
model_0_5_lower = getModelOutput_matlab([kPB_lower, F_lower], vector_time_points, 0.5, Vol_P, Vol_B, k2);

model_1_0_upper = getModelOutput_matlab([kPB_upper, F_upper], vector_time_points, 1.0, Vol_P, Vol_B, k2);
model_1_0_lower = getModelOutput_matlab([kPB_lower, F_lower], vector_time_points, 1.0, Vol_P, Vol_B, k2);

%% ---------------------- Figure ----------------------
figure('Color', 'w'); hold on;

hCI05 = fill_ci_band(vector_time_points, model_0_5_lower, model_0_5_upper, col_model_05, ci_alpha);
hCI10 = fill_ci_band(vector_time_points, model_1_0_lower, model_1_0_upper, col_model_10, ci_alpha);

hM05 = plot(vector_time_points, model_0_5, '-', ...
    'Color', col_model_05, 'LineWidth', 1.8);

hM10 = plot(vector_time_points, model_1_0, '-', ...
    'Color', col_model_10, 'LineWidth', 1.8);

hD05 = plot_asym_errorbars(time_points, dose_0_5, err05_lower, err05_upper, col_data_05, 6);
hD10 = plot_asym_errorbars(time_points, dose_1_0, err10_lower, err10_upper, col_data_10, 6);

xlabel('Time (hours)');
ylabel('LPS Concentration (ng/mL)');
title('Pharmacokinetic Subsystem Fitting (Blood)');

grid off;
xticks(0:2:24);
xlim([-1, 25]);

set(gca, 'LineWidth', 1.5, 'FontWeight', 'bold', 'TickDir', 'out');
box off;

legend([hD05, hD10, hM05, hM10, hCI05, hCI10], ...
       {'Observed data 0.5 mg/kg', ...
        'Observed data 1.0 mg/kg', ...
        'Model 0.5 mg/kg', ...
        'Model 1.0 mg/kg', ...
        'Model 95% CI 0.5 mg/kg', ...
        'Model 95% CI 1.0 mg/kg'}, ...
       'FontSize', 8, 'Location', 'best');

%% ====================== Local functions ======================

function y = FObjective_data_space(params, time_points, Vol_P, Vol_B, k2)

    B05 = getModelOutput_matlab(params, time_points, 0.5, Vol_P, Vol_B, k2);
    B10 = getModelOutput_matlab(params, time_points, 1.0, Vol_P, Vol_B, k2);

    y = [B05; B10];

end

function B_sol = getModelOutput_matlab(params, t_eval, dose, Vol_P, Vol_B, k2)

    kPB = params(1);
    F   = params(2);

    P0 = dose * 248000;
    B0 = 0;

    y0 = [P0; B0];

    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    [t, y] = ode45(@(t,y) systemODE_matlab(t, y, kPB, F, Vol_P, Vol_B, k2), ...
                   [0 24], y0, opts);

    B_sol = interp1(t, y(:,2), t_eval, 'pchip');

    % Prevent accidental log(0) issues during fitting.
    B_sol = max(B_sol, realmin);

end

function dydt = systemODE_matlab(~, Y, kPB, F, Vol_P, Vol_B, k2)

    P = Y(1);
    B = Y(2);

    dPdt = -kPB * P;
    dBdt = F * kPB * (Vol_P / Vol_B) * P - k2 * B;

    dydt = [dPdt; dBdt];

end

function h = plot_asym_errorbars(x, y, elow, ehigh, col, markersize)

    h = plot(x, y, 'o', ...
        'MarkerSize', markersize, ...
        'LineWidth', 1.0, ...
        'MarkerEdgeColor', col, ...
        'MarkerFaceColor', 'w', ...
        'Color', col);

    xr = max(x) - min(x);
    cap = max(0.02 * xr, 0.15);

    for i = 1:numel(x)
        line([x(i) x(i)], [y(i)-elow(i), y(i)+ehigh(i)], ...
            'Color', col, 'LineWidth', 1.0);

        line([x(i)-cap x(i)+cap], [y(i)-elow(i), y(i)-elow(i)], ...
            'Color', col, 'LineWidth', 1.0);

        line([x(i)-cap x(i)+cap], [y(i)+ehigh(i), y(i)+ehigh(i)], ...
            'Color', col, 'LineWidth', 1.0);
    end

end

function h = fill_ci_band(x, ylow, yhigh, col, alphaVal)

    h = fill([x; flipud(x)], [ylow; flipud(yhigh)], col, ...
        'LineStyle', 'none');

    set(h, 'FaceAlpha', alphaVal);
    uistack(h, 'bottom');

end

function J = finiteDiffJacobian(fun, p0)

    f0 = fun(p0);
    m = numel(f0);
    n = numel(p0);

    J = zeros(m, n);

    h = 1e-6 * max(1, abs(p0(:)));

    for j = 1:n
        pj1 = p0;
        pj2 = p0;

        pj1(j) = pj1(j) + h(j);
        pj2(j) = pj2(j) - h(j);

        J(:,j) = (fun(pj1) - fun(pj2)) / (2 * h(j));
    end

end