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
% CDF_hat:              grid for the aproximated CDF

N = 2^Mfft; % Needed for optimal du computation
[du, a] = extraParamsComputation(model, activity, params, N, dt, STRETCH);

if params(1)>0 && params(1) < 1 && strcmp(model, 'NTS-OU')
    du = du*params(1);
end

% Compute the parameters for the FFT discretization.
u_N = 0.5 * du * (N - 1); u_1 = -u_N;
dx = 2 * pi / (N * du); x_1 = -dx * (N - 1) / 2; x_N = -x_1;

Ra = a < 0;

% Assign values to the struct for the FFT function.
numericalParams.M = Mfft;
numericalParams.u1 = u_1;
numericalParams.uN = u_N;
numericalParams.du = du;
numericalParams.x1 = x_1;
numericalParams.xN = x_N;
numericalParams.dx = dx;

% Filter the xgrid to avoid errors when selecting the largest subset.
xgrid = x_1:dx:x_N;
xgrid = xgrid(xgrid >= -20 & xgrid <= 20);

% FFT to retrieve the CDF on the xgrid.
% First assign the function for the FFT.
f = @(u) exp(LogCharFunc(STRETCH*(u + 1i .* a), dt, params, model, activity))./ (1i .* (u + 1i .* a));

I_fft = computeFFT(f, xgrid, numericalParams);
RawCDF = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft;

% Remove NaN values.
valid_idxs = ~isnan(RawCDF);
xgrid = xgrid(valid_idxs);
RawCDF = RawCDF(valid_idxs);

% Restrict the xgrid by taking the largest set where the raw CDF is
% monotone increasing.
[xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid);
[CDF_hat, idxs_unique] = unique(CDF_hat);
xgrid_hat = xgrid_hat(idxs_unique);
end