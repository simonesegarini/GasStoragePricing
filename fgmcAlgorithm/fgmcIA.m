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

% Plot of the CDF obtained by inverting the CF and of the final plot with
% all the infos needed for debugging, uncomment if probelms arise.
%
% alpha = params(1);
%
% figure;
% plot(xgrid, RawCDF, '-k');
% title(['Raw CDF with alpha = ', num2str(alpha)])
% 
% disp(['Caso: ', model, ', alpha = ', num2str(alpha), ', du = ', num2str(du)])
% disp(['Extrapolated points: ', num2str(length(U)-sum(idxs_within))])
% disp(' ')
% 
% figure;
% plot(xgrid_hat, CDF_hat, '--r')
% title(['Plot with alpha = ', num2str(alpha)])
% hold on;
% plot(increments(1:100), U(1:100), 'ok');
% first_last = [xgrid_hat(1), CDF_hat(1); xgrid_hat(end), CDF_hat(end)];
% plot(first_last(:, 1), first_last(:, 2), '*g')
% plot(incrOver, UOver, 'db')
% plot(incrUnder, UUnder, 'db')
% legend('Interpolated Inverted CDF', 'Interpolated values', 'Limits', 'Exponential extrapolation Over', 'Exponential extrapolation Under')
end