function increments = fgmcIA(U, params, Mfft, dt, model, activity, toll)
% FGMC method for Infinite Activity processes
%
% INPUT:
% U:                    matrix with the simulated Uniform RV
% params:               vector of paramteres for the specified model
% M:                    parameter for FFT
% dt:                   time step
% model:                model selected
% activity:             model activity, needed for the FA case
% toll:                 tollerance for CDF selection
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy

N = 2^Mfft; % Needed for optimal du computation
[du, a] = extraParamsComputation(model, activity, params, N, dt);

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

% Filter the xgrid to avoid errors when selecting the largest subset
xgrid = x_1:dx:x_N;
xgrid = xgrid(xgrid >= -20 & xgrid <= 20);

% FFT to retrieve the CDF on the xgrid.
% First assign the function for the FFT.
f = @(u) exp(LogCharFunc(u + 1i .* a, dt, params, model, activity)) ./ (1i .* (u + 1i .* a));

I_fft = computeFFT(f, xgrid, numericalParams);
RawCDF = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft;

% Remove NaN values.
valid_idxs = ~isnan(RawCDF);
xgrid = xgrid(valid_idxs);
RawCDF = RawCDF(valid_idxs);

% Restrict the xgrid by taking the largest set where the raw CDF is
% monotone increasing.
[xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid, toll);
[CDF_hat, idxs_unique] = unique(CDF_hat);
xgrid_hat = xgrid_hat(idxs_unique);

% Use spline interpolation within the range and exponential extrapolation
% outside.
increments = zeros(size(U));
idxs_within = U >= CDF_hat(1) & U <= CDF_hat(end);
increments(idxs_within) = interp1(CDF_hat, xgrid_hat, U(idxs_within), 'spline');

% Exponential extrapolation, dx as U_n+1 = 1-exp(-b*X_n+1), sx as U_n-1 = exp(b*X_n-1)
% Search for the values that have to be extrapolated.
UOver = U(U>CDF_hat(end));
UUnder = U(U<CDF_hat(1));

bOver = -1/xgrid_hat(end) * log(1-CDF_hat(end));
incrOver = -1/bOver*log(1-UOver);

bUnder = -1/xgrid_hat(1)*log(CDF_hat(1));
incrUnder = -1/bUnder*log(UUnder);

% Assign the values extrapolated.
increments(U > CDF_hat(end)) = incrOver;
increments(U < CDF_hat(1)) = incrUnder;

end