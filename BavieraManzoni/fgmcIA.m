function increments = fgmcIA(U, params, M, dt, model, activity)
% FGMC method for Infinite Activity processes
%
% INPUT:
% U:                    matrix with the simulated Uniform RV
% params:               vector of paramteres for the specified model
% M:                    parameter for FFT
% dt:                   time step
% model:                model selected
% activity:             model activity, needed for the FA case
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy

% Proceed to assign parameters, compute analyticity strip (p_n, p_p) and
% compute the optimal step du for the FFT code.
% For the finite activity processes du was found by checking the one that
% gave empirical cumulants closer to the theoretical ones.
%
% Finite activity processes have power decay, infinite activity have
% exponential decay. Special cases are commented.

alpha = params(1); b = params(2);
N = 2^M; % Needed for optimal du computation

switch model
    case {'OU-NTS', 'NTS-OU'}
        sigma = params(3);
        k = params(4); theta = params(5);
    
        A_as = sqrt(theta^2 + (2 * sigma^2 * (1 - alpha)) / k);
        p_n = (theta - A_as) / (sigma^2);
        p_p = (theta + A_as) / (sigma^2);
        
        switch activity
            case 'Finite'
                a = abs(0.25 * max(-p_n, p_p));
                du = 0.0971; % Best for NTS
            case 'Infinite'
                if alpha == 0 % OU-VG case (Variance Gamma, NTS), power decay.
                    a = 0.25 * max(-p_n, p_p);
                    du = 0.0971;
                else
                    a = 0.5 * max(-p_n, p_p);
                    omega = 2 * alpha;
                    if strcmp(model, 'OU-NTS') % OU-NTS GENERAL CASE
                        l = 1 / alpha * ((1 - alpha) / k)^(1 - alpha) * (sigma^2 / 2)^alpha * (1 - exp(-2 * alpha * b * dt)) / (2 * alpha * b);
                    else % NTS-OU GENERAL CASE
                        l = 1 / alpha * ((1 - alpha) / k)^(1 - alpha) * (sigma^2 / 2)^alpha * (1 - exp(-2 * alpha * b * dt));
                    end
                    du = (2 * pi * abs(a) / (l * N^omega))^(1 / (omega + 1));
                end
        end
        
    case {'OU-TS', 'TS-OU'}
        beta_p = params(3); beta_n = params(4); c_p = params(5); c_n = params(6);
    
        p_n = -beta_p;
        p_p = beta_n;

        switch activity
            case 'Finite'
                a = abs(0.25 * max(-p_n, p_p));
                du = 0.4; % Best for TS
            case 'Infinite'
                if alpha == 0 % OU-Gamma case (TS), power decay.
                    a = 0.25 * max(-p_n, p_p);
                    du = 0.0971;
                else
                    a = 0.5 * max(-p_n, p_p);
                    if strcmp(model, 'OU-TS')
                        if alpha == 1 % Other special case for OU-TS, exp decay.
                            omega = 1;
                            l = (c_p + c_n) * pi / 2 * (1 - exp(-b * dt)) / b;
                        else % GENERAL OU-TS CASE.
                            omega = alpha;
                            scale = 1;
                            l = -(c_p + c_n) * gamma(-alpha) * cos(alpha * pi / 2) * (1 - exp(-alpha * b * dt)) / (alpha * b) * scale^alpha;
                        end
                    else % TS-OU
                        if alpha == 1 % Other special case for TS-OU, exp decay.
                            omega = 1;
                            l = (c_p + c_n) * pi / 2 * (1 - exp(-b * dt));
                        else % GENERAL TS-OU CASE.
                            omega = alpha;
                            scale = 1;
                            l = -(c_p + c_n) * gamma(-alpha) * cos(alpha * pi / 2) * (1 - exp(-alpha * b * dt)) * scale^alpha;
                        end
                    end
                    du = (2 * pi * abs(a) / (l * N^omega))^(1 / (omega + 1));
                end
        end
end

% Compute the parameters for the FFT discretization.
u_N = 0.5 * du * (N - 1); u_1 = -u_N;
dx = 2 * pi / (N * du); x_1 = -dx * (N - 1) / 2; x_N = -x_1;

Ra = a < 0;

disp(['Caso: ', model, ', alpha = ', num2str(alpha), ', du = ', num2str(du)])

% Assign values to the struct for the FFT function.
numericalParams.M = M;
numericalParams.u1 = u_1;
numericalParams.uN = u_N;
numericalParams.du = du;
numericalParams.x1 = x_1;
numericalParams.xN = x_N;
numericalParams.dx = dx;

xgrid = -20:dx:20;

% FFT to retrieve the CDF on the xgrid.
% First assign the function for the FFT.
f = @(u) exp(LogCharFunc(u + 1i .* a, dt, params, model, activity)) ./ (1i .* (u + 1i .* a));

I_fft = computeFFT(f, xgrid, numericalParams);
RawCDF = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft;

% Remove NaN values.
valid_idxs = ~isnan(RawCDF);
xgrid = xgrid(valid_idxs);
RawCDF = RawCDF(valid_idxs);

% Plot of the CDF obtained by inverting the CF, uncomment if needed.
% figure;
% plot(xgrid, RawCDF, '-k');
% title(['Raw CDF with alpha = ', num2str(alpha)])

% Restrict the xgrid by taking the largest set where the raw CDF is
% monotone increasing.
[xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid, 1e-9);
[CDF_hat, idxs_unique] = unique(CDF_hat);
xgrid_hat = xgrid_hat(idxs_unique);


% Use spline interpolation within the range and exponential extrapolation
% outside.
increments = zeros(size(U));
idxs_within = U >= CDF_hat(1) & U <= CDF_hat(end);
increments(idxs_within) = interp1(CDF_hat, xgrid_hat, U(idxs_within), 'spline');
disp(['Extrapolated points: ', num2str(length(U)-sum(idxs_within))])
disp(' ')

% Exponential extrapolation
% dx as U_n+1 = 1-a*exp(-b*X_n+1)
% sx as U_n-1 = a*exp(b*X_n-1)

UOver = U(U>CDF_hat(end));
UUnder = U(U<CDF_hat(1));

bOver = log((1-CDF_hat(end-1))/(1-CDF_hat(end)))/(xgrid_hat(end)-xgrid_hat(end-1));
aOver = (1-CDF_hat(end))*exp(bOver*xgrid_hat(end));

bUnder = log(CDF_hat(2)/CDF_hat(1))/(xgrid_hat(2)-xgrid_hat(1));
aUnder = exp(-bUnder*xgrid_hat(1))*CDF_hat(1);

incrOver = -1/bOver .* log((1-UOver)./aOver);
incrUnder = 1/bUnder .* log(UUnder./aUnder);

% aUnder = CDF_hat(1);
% bUnder = log(CDF_hat(2)/CDF_hat(1))/(xgrid_hat(2)-xgrid_hat(1));
% 
% aOver = CDF_hat(end);
% bOver = log((1-CDF_hat(end-1))/(1-CDF_hat(end))) / (xgrid_hat(end-1)-xgrid_hat(end));
% 
% incrUnder = log(UUnder/aUnder)/bUnder + xgrid_hat(1);
% incrOver = log(-(UOver-1)/aOver)/bOver+xgrid_hat(end);

% Assign the values extrapolated.
increments(U > CDF_hat(end)) = incrOver;
increments(U < CDF_hat(1)) = incrUnder;

% Final plot, uncomment if needed.
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