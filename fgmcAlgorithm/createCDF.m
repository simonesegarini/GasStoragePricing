function [xgrid_hat, CDF_hat] = createCDF(params, Mfft, dt, model, activity, STRETCH)
% Function that creates the approximated CDF of the increments.
%
% INPUT:
% params:               vector with the parameters of the model
% Mfft:                 FFT hyperparameter
% dt:                   time horizon
% model:                char for model selection
% activity:             char for model selection
% STRETCH:              parameter to stretch the CDF and handle low alphas
%
% OUTPUT:
% xgrid_hat:            grid for the x values
% CDF_hat:              grid for the approximated CDF

N = 2^Mfft; % Number of FFT points
[du, a] = extraParamsComputation(model, activity, params, N, dt, STRETCH);

if params(1)>0 && params(1) < 1
    du = du*params(1);
end

% Compute the parameters for the FFT discretization.
dx = 2 * pi / (N * du); 
x_1 = -dx * (N - 1) / 2; x_N = -x_1;
u_N = 0.5 * du * (N - 1); u_1 = -u_N;

Ra = a < 0;
disp(['Step du = ', num2str(du)])
disp(['Step dx = ', num2str(dx)])
disp(['Params a = ', num2str(a)])
disp(' ')

% Assign values to the struct for the FFT function.
numericalParams.M = Mfft;
numericalParams.u1 = u_1;
numericalParams.uN = u_N;
numericalParams.du = du;
numericalParams.x1 = x_1;
numericalParams.xN = x_N;
numericalParams.dx = dx;

% Define the xgrid and apply scaling
xgrid = x_1:dx:x_N;
xgrid = xgrid(xgrid >= -20 & xgrid <= 20);

% Define the characteristic function with the scaling
f = @(u) exp((LogCharFunc(STRETCH*(u + 1i * a), dt, params, model, activity))) ./ (1i *(u + 1i * a));

% Compute FFT
I_fft = computeFFT(f, xgrid, numericalParams);
RawCDF = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft;

% Remove NaN values
valid_idxs = ~isnan(RawCDF);
xgrid = xgrid(valid_idxs);
RawCDF = RawCDF(valid_idxs);

% Restrict the xgrid by taking the largest set where the raw CDF is monotone increasing
[xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid);
[CDF_hat, idxs_unique] = unique(CDF_hat);
xgrid_hat = xgrid_hat(idxs_unique);

% Plot for debugging
alpha = params(1);
figure;
plot(xgrid, RawCDF, '-k');
title(['Plot with alpha = ', num2str(alpha)])
xlim([-1, 1])
hold on
plot(xgrid_hat, CDF_hat, '--r')
first_last = [xgrid_hat(1), CDF_hat(1); xgrid_hat(end), CDF_hat(end)];
plot(first_last(:, 1), first_last(:, 2), '*g')
ylim([0, 1])
end
