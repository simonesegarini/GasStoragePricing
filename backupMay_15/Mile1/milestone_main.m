%% MILESTONE 1
clear all, close all, clc

%% SIMULATION PARAMETERS
T = 1/12;
M = 1;
N = 1e7;

%% SPOT PARAMETERS
X0 = 0;
F0t = 1; %flat forward curve

%% OU-NTS
alpha = 0.4; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;
tic
[Xt, ft]= fgmc(X0, [alpha, b, sigma, k, theta], N, M, T, 'NTS');
time = toc;
[TCumulantsNTS, ECumulantsNTS] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'NTS')

S = F0t.*exp(Xt + ft);

%% OU-NTS VARYING ALPHA 
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

TCumulantsNTSalpha = zeros(numel(alphas), 4);
ECumulantsNTSalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    [Xt, ~]= fgmc(X0, [alphas(i), b, sigma, k, theta], N, M, T, 'NTS');
    [TCumulantsNTSalpha(i,:), ECumulantsNTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'NTS');
end
%% OU-TS, TO FIX
alpha = 1.6; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

[Xt, ft] = fgmc(X0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'TS');
[TCumulantsTS, ECumulantsTS] = computeCumulants(Xt, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS')

S = F0t.*exp(Xt + ft);

%% OU-TS VARYING ALPHA,TO FIX
alphas = [1.6, 1.2, 0.8, 0.4, -1.0, -2.0];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSalpha = zeros(numel(alphas), 4);
ECumulantsTSalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    [Xt, ~]= fgmc(X0, [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'TS');
    [TCumulantsTSalpha(i,:), ECumulantsTSalpha(i,:)] = computeCumulants(Xt(:, end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS');
end