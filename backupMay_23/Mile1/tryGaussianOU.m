clear all, close all, clc

%% SIMULATION PARAMETERS
T = 1/12;
M = 1;
N = 1e7;

%% SPOT PARAMETERS
X0 = 0;
F0t = 1; %flat forward curve

%%
rng(1) % set seed for reproducibility
dt = T/M;
U = rand(N, M);
X = zeros(N,M+1); X(:,1) = X0;
% f = zeros(1,M+1); f(:,1) = -LogCharFunc(-1i, 0, params, model);
Mfft = 16;
alpha = 0.8;
b = 0.2162;

integrand = @(z) exp(-2*b*(dt-z));
z_values = linspace(0, dt, 1000);
sigma = sqrt(trapz(z_values, integrand(z_values)));

%%
charFunc = @(u) exp(-0.5.*u.^2.*sigma.^2);

N=2^Mfft;
u_1 = -500; u_N=-u_1; du=(u_N-u_1)/(N-1); 
dx=2*pi/(N*du); x_1=-dx*(N-1)/2; x_N=-x_1;

numericalParams.M=Mfft;
numericalParams.u1=u_1;
numericalParams.uN=u_N;
numericalParams.du=du;
numericalParams.x1=x_1;
numericalParams.xN=x_N;
numericalParams.dx=dx;

xgrid = -5:dx:5;
%% Iteration in discretized time to update the simulations

f = @(u) charFunc(u)./(1i.*u);
tic
for j = 1:M
    
    [I_fft, ~] = computeIntegral(f, xgrid, numericalParams);
    RawCDF = 0.5 - I_fft./(2*pi);

    figure;
    plot(xgrid, RawCDF, '-k');
    title(['Raw CDF with alpha = ', num2str(alpha)])

    [xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid);

    % Remove duplicate CDF_hat values and corresponding xgrid values
    [CDF_hat_unique, idx] = unique(CDF_hat);
    xgrid_unique = xgrid_hat(idx);
    
    % Invert the CDF using a spline interpolation
    inverse_interpolator = @(y) interp1(CDF_hat_unique, xgrid_unique, y, 'spline');
    increments = inverse_interpolator(U(:,j));

    figure;
    plot(xgrid, RawCDF, '-b');
    title(['Plot with alpha =  ', num2str(alpha)])
    hold on;
    plot(xgrid_hat, CDF_hat, '--r')
    plot(increments(1:100), U(1:100,j), 'ok');
    legend('Inverted CDF', 'Interpolated Inverted CDF', 'Interpolated values')


    X(:, j+1) = exp(-b.*dt).*X(:, j) + increments;
%     f(:, j+1) = - LogCharFunc(-1i, dt*j, params, model); 
end
time = toc;

%% plot just to check
figure;
hold on

for i=1:20
    plot(exp(X(i,:)));
end
