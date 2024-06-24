clear all, close all, clc

% Define your data points
xgrid = [1, 2, 3, 4, 5];
ygrid = [0.5, 0.7, 0.85, 0.93, 0.97];

% Separate data into two parts based on a threshold (e.g., median or a chosen value)
threshold = 0.85;
above_indices = ygrid > threshold;
below_indices = ygrid <= threshold;

x_above = xgrid(above_indices);
y_above = ygrid(above_indices);
x_below = xgrid(below_indices);
y_below = ygrid(below_indices);

% Define the exponential functions
expFuncAbove = @(params, x) 1 - params(1) * exp(-params(2) * x);
expFuncBelow = @(params, x) params(1) * exp(params(2) * x);

% Fit the curve to the data above the threshold
params0_above = [1, 1]; % Initial guess for parameters [a1, b1]
params_above = lsqcurvefit(expFuncAbove, params0_above, x_above, y_above);

% Fit the curve to the data below the threshold
params0_below = [1, 1]; % Initial guess for parameters [a2, b2]
params_below = lsqcurvefit(expFuncBelow, params0_below, x_below, y_below);

% Extract fitted parameters
a1 = params_above(1);
b1 = params_above(2);
a2 = params_below(1);
b2 = params_below(2);

% Values to extrapolate above and below the threshold
u_values_above = [0.98, 0.99, 0.995];
u_values_below = [0.4, 0.3, 0.2];

% Extrapolate x for these u values above the threshold
x_extrapolated_above = -log((1 - u_values_above) / a1) / b1;

% Extrapolate x for these u values below the threshold
x_extrapolated_below = log(u_values_below / a2) / b2;

% Display the results
disp(['Fitted parameters above threshold: a1 = ', num2str(a1), ', b1 = ', num2str(b1)]);
disp(['Extrapolated x values for u above threshold: ', num2str(x_extrapolated_above)]);

disp(['Fitted parameters below threshold: a2 = ', num2str(a2), ', b2 = ', num2str(b2)]);
disp(['Extrapolated x values for u below threshold: ', num2str(x_extrapolated_below)]);

% Plotting the data and the extrapolation
x_fit = linspace(min(xgrid), max(xgrid) + 2, 100);
y_fit_above = expFuncAbove(params_above, x_fit);
y_fit_below = expFuncBelow(params_below, x_fit);

figure;
hold on;
scatter(xgrid, ygrid, 'ro', 'DisplayName', 'Data');
plot(x_fit, y_fit_above, 'b-', 'DisplayName', 'Fitted exponential curve (Above)');
plot(x_fit, y_fit_below, 'g-', 'DisplayName', 'Fitted exponential curve (Below)');
scatter(x_extrapolated_above, u_values_above, 'bo', 'DisplayName', 'Extrapolated points (Above)');
scatter(x_extrapolated_below, u_values_below, 'go', 'DisplayName', 'Extrapolated points (Below)');
xlabel('x');
ylabel('y');
legend('show');
hold off;
