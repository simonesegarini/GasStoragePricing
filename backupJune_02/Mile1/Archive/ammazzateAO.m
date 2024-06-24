clear all, close all, clc

warningID = 'MATLAB:legend:IgnoringExtraEntries';
warning('off', warningID);

% SIMULATION PARAMETERS
T = 1/12;
M = 1;
N = 1e7;

% SPOT PARAMETERS
X0 = 0;

rng(2) % set seed for reproducibility
dt = T/M;
U = rand(N, M);
X_1 = zeros(N,M+1); X_1(:,1) = X0;
X_2 = zeros(N,M+1); X_2(:,1) = X0;

%% OU-NTS
alpha = 0.8; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

A_as = sqrt(theta^2 + (2*sigma^2*(1-alpha)) / k);
p_n = (theta - A_as) / (sigma.^2);
p_p = (theta + A_as) / (sigma.^2);
a = 0.5 * max(-p_n, p_p);

N = 2^16;
omega = 2 * alpha;
l = 1 / alpha * ((1 - alpha) / k)^(1 - alpha) * (sigma^2 / 2)^alpha * (1 - exp(-2 * alpha * b * dt)) / (2 * alpha * b);
Ra = a < 0;


du = (2 * pi * abs(a) / (l * N^omega))^(1 / (omega + 1));
u_1 = -(N - 1) / 2 * du; u_N = -u_1;
dx = 2 * pi / (N * du);
x_1 = -dx * (N - 1) / 2; x_N = -x_1;

numericalParams.M = 16;
numericalParams.u1 = u_1;
numericalParams.uN = u_N;
numericalParams.du = du;
numericalParams.x1 = x_1;
numericalParams.xN = x_N;
numericalParams.dx = dx;

xgrid = -5:dx:5;

[I_fft_1, I_fft_2] = computeIntegral(xgrid, numericalParams, [alpha, b, sigma, k, theta], a, dt);
RawCDF_1 = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft_1;
RawCDF_2 = Ra - exp(a .* xgrid) ./ (2 * pi) .* I_fft_2;

figure;
plot(xgrid, RawCDF_1, '-k');
figure;
plot(xgrid, RawCDF_2, '-k');

% Restrict the xgrid by taking the largest set where the raw CDF is monotone increasing
[xgrid_hat_1, CDF_hat_1] = approxCDF(RawCDF_1, xgrid);
[xgrid_hat_2, CDF_hat_2] = approxCDF(RawCDF_2, xgrid);

% Use spline interpolation within the range and exponential extrapolation outside
increments_1 = zeros(size(U));
idxs_within_1 = U >= CDF_hat_1(1) & U <= CDF_hat_1(end);
increments_1(idxs_within_1) = interp1(CDF_hat_1, xgrid_hat_1, U(idxs_within_1), 'spline');

UOver_1 = U(U>CDF_hat_1(end));
UUnder_1 = U(U<CDF_hat_1(1));

bOver_1 = log((1-CDF_hat_1(end-1))/(1-CDF_hat_1(end)))/(xgrid_hat_1(end)-xgrid_hat_1(end-1));
aOver_1 = (1-CDF_hat_1(end))*exp(bOver_1*xgrid_hat_1(end));

bUnder_1 = log(CDF_hat_1(2)/CDF_hat_1(1))/(xgrid_hat_1(2)-xgrid_hat_1(1));
aUnder_1 = exp(-bUnder_1*xgrid_hat_1(1))*CDF_hat_1(1);

incrOver_1 = -1/bOver_1 .* log((1-UOver_1)./aOver_1);
incrUnder_1 = 1/bUnder_1 .* log(UUnder_1./aUnder_1);

increments_1(U > CDF_hat_1(end)) = incrOver_1;
increments_1(U < CDF_hat_1(1)) = incrUnder_1;

increments_2 = zeros(size(U));
idxs_within_2 = U >= CDF_hat_2(1) & U <= CDF_hat_2(end);
increments_2(idxs_within_2) = interp1(CDF_hat_2, xgrid_hat_2, U(idxs_within_2), 'spline');

UOver_2 = U(U>CDF_hat_2(end));
UUnder_2 = U(U<CDF_hat_2(1));

bOver_2 = log((1-CDF_hat_2(end-1))/(1-CDF_hat_2(end)))/(xgrid_hat_2(end)-xgrid_hat_2(end-1));
aOver_2 = (1-CDF_hat_2(end))*exp(bOver_2*xgrid_hat_2(end));

bUnder_2 = log(CDF_hat_2(2)/CDF_hat_2(1))/(xgrid_hat_2(2)-xgrid_hat_2(1));
aUnder_2 = exp(-bUnder_2*xgrid_hat_2(1))*CDF_hat_2(1);

incrOver_2 = -1/bOver_2 .* log((1-UOver_2)./aOver_2);
incrUnder_2 = 1/bUnder_2 .* log(UUnder_2./aUnder_2);

increments_2(U > CDF_hat_2(end)) = incrOver_2;
increments_2(U < CDF_hat_2(1)) = incrUnder_2;

figure;
plot(xgrid_hat_1, CDF_hat_1, '--r')
title(['Plot with alpha = ', num2str(alpha)])
hold on;
plot(increments_1(1:100), U(1:100), 'ok');
first_last = [xgrid_hat_1(1), CDF_hat_1(1); xgrid_hat_1(end), CDF_hat_1(end)];
plot(first_last(:, 1), first_last(:, 2), '*g')
plot(incrOver_1, UOver_1, 'db')
plot(incrUnder_1, UUnder_1, 'db')
legend('Interpolated Inverted CDF', 'Interpolated values', 'Limits', 'Exponential extrapolation Over', 'Exponential extrapolation Under')

figure;
plot(xgrid_hat_2, CDF_hat_2, '--r')
title(['Plot with alpha = ', num2str(alpha)])
hold on;
plot(increments_2(1:100), U(1:100), 'ok');
first_last = [xgrid_hat_2(1), CDF_hat_2(1); xgrid_hat_2(end), CDF_hat_2(end)];
plot(first_last(:, 1), first_last(:, 2), '*g')
plot(incrOver_2, UOver_2, 'db')
plot(incrUnder_2, UUnder_2, 'db')
legend('Interpolated Inverted CDF', 'Interpolated values', 'Limits', 'Exponential extrapolation Over', 'Exponential extrapolation Under')

X_1(:, 2) = exp(-b.*dt).*X_1(:, 1) + increments_1;
X_2(:, 2) = exp(-b.*dt).*X_2(:, 1) + increments_2;

[TCumulantsOUNTS_1, ECumulantsOUNTS_1] = computeCumulants(X_1(:,end), [alpha, b, sigma, k, theta], T, 'OU-NTS')
[TCumulantsOUNTS_2, ECumulantsOUNTS_2] = computeCumulants(X_2(:,end), [alpha, b, sigma, k, theta], T, 'OU-NTS')
