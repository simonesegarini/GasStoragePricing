function increments = fgmcIA(xgrid_hat, CDF_hat, U)
% FGMC method for Infinite Activity processes.
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

% disp(['Extrapolated points: ', num2str(length(U)-sum(idxs_within))])
% disp(' ')
end