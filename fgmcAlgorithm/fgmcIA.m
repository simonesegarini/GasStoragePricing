function increments = fgmcIA(xgrid_hat, CDF_hat, U)
% FGMC method for Infinite Activity processes.
%
% INPUT:
% xgrid_hat:            xgrid of reconstructed cdf
% CDF_hat               ygrid of reconstructed cdf
% U:                    matrix with the simulated Uniform RV
% STRETCH:              scaling factor used in CDF computation
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy

% Use spline interpolation within the range and exponential extrapolation
% outside.
increments = zeros(size(U));
idxs_within = U >= CDF_hat(1) & U <= CDF_hat(end);
if any(idxs_within)
    increments(idxs_within) = interp1(CDF_hat, xgrid_hat, U(idxs_within), 'spline');
end

% Exponential extrapolation, adjust based on scaling
UOver = U(U > CDF_hat(end));
if ~isempty(UOver)
    bOver = -1 / xgrid_hat(end) * log(1 - CDF_hat(end));
    incrOver = -1 / bOver * log(1 - UOver);
    increments(U > CDF_hat(end)) = incrOver;
end

UUnder = U(U < CDF_hat(1));
if ~isempty(UUnder)
    bUnder = -1 / xgrid_hat(1) * log(CDF_hat(1));
    incrUnder = -1 / bUnder * log(UUnder);
    increments(U < CDF_hat(1)) = incrUnder;
end

% disp(['Extrapolated points: ', num2str(length(U) - sum(idxs_within))])
% disp(' ')
end
