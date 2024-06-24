%% MILESTONE 1
clear all, close all, clc

warningID = 'MATLAB:legend:IgnoringExtraEntries';
warning('off', warningID);

%% SIMULATION PARAMETERS
T = 1/12;
M = 1;
N = 1e7;

%% SPOT PARAMETERS
X0 = 0;

%% OU-NTS
alpha = 0.2; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;
tic
Xt = fgmc(X0, [alpha, b, sigma, k, theta], N, M, T, 'OU-NTS');
time = toc;
[TCumulantsOUNTS, ECumulantsOUNTS] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'OU-NTS')

%% OU-NTS VARYING ALPHA 
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

TCumulantsOUNTSalpha = zeros(numel(alphas), 4);
ECumulantsOUNTSalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    Xt = fgmc(X0, [alphas(i), b, sigma, k, theta], N, M, T, 'OU-NTS');
    [TCumulantsOUNTSalpha(i,:), ECumulantsOUNTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'OU-NTS');
end

%% NTS-OU
alpha = 0.2; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;
tic
Xt = fgmc(X0, [alpha, b, sigma, k, theta], N, M, T, 'NTS-OU');
time = toc;
[TCumulantsNTSOU, ECumulantsNTSOU] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'NTS-OU')

%% NTS-OU VARYING ALPHA 
alphas = [0.8, 0.6, 0.4, 0.2];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

TCumulantsNTSOUalpha = zeros(numel(alphas), 4);
ECumulantsNTSOUalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    Xt = fgmc(X0, [alphas(i), b, sigma, k, theta], N, M, T, 'NTS-OU');
    [TCumulantsNTSOUalpha(i,:), ECumulantsNTSOUalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'NTS-OU');
end


%% OU-TS, TO FIX
alpha = 0.4; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

Xt = fgmc(X0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'OU-TS');
[TCumulantsOUTS, ECumulantsOUTS] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS')

%% OU-TS VARYING ALPHA,TO FIX
alphas = [1.6, 1.2, 0.8, 0.4, -1.0, -2.0];
%alphas = [1.6, 1.2, 0.8, 0.4];
%alphas = [-1.0, -2.0];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsOUTSalpha = zeros(numel(alphas), 4);
ECumulantsOUTSalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    Xt= fgmc(X0, [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'OU-TS');
    [TCumulantsOUTSalpha(i,:), ECumulantsOUTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS');
end

%% TS-OU, TO FIX
alpha = 1.6; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

Xt = fgmc(X0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'TS-OU');
[TCumulantsTS0U, ECumulantsTS0U] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU')

%% TS-OU VARYING ALPHA,TO FIX
alphas = [1.6, 1.2, 0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSOUalpha = zeros(numel(alphas), 4);
ECumulantsTSOUalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    Xt= fgmc(X0, [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 'TS-OU');
    [TCumulantsTSOUalpha(i,:), ECumulantsTSOUalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU');
end