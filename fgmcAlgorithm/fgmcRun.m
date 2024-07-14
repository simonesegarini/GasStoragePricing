%% fgmc - IMPLEMENTATION OF BAVIERA MANZONI 2024 PAPER
clear all, close all, clc

% Turn off the warning for extra legend entries, just to avoid messages
% when using the legend function if there are no extrapolated values to
% plot. Useless if the plor in fgmcIA are commented.
warningID = 'MATLAB:legend:IgnoringExtraEntries';
warning('off', warningID);

%% SIMULATION PARAMETERS
T = 1; % Time horizon
M = 1; % Number of steps
N = 1e7; % Number of simulations
seed = 2; % Seed for the random uniform sampling

%% Processes listed separately to run with the alphas of the paper all together

%% OU-NTS VARYING ALPHA, SYMMETRIC
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

TCumulantsOUNTSalphaSymm = zeros(numel(alphas), 4);
ECumulantsOUNTSalphaSymm = zeros(numel(alphas), 4);
timeOUNTSalphaSymm = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt = fgmc('OU-NTS', [alphas(i), b, sigma, k, theta], N, M, T, 16, seed, 1e-8);
    [TCumulantsOUNTSalphaSymm(i,:), ECumulantsOUNTSalphaSymm(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'OU-NTS');
end

%% OU-NTS VARYING ALPHA, ASYMMETRIC
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

TCumulantsOUNTSalphaAsymm = zeros(numel(alphas), 4);
ECumulantsOUNTSalphaAsymm = zeros(numel(alphas), 4);
timeOUNTSalphaAsymm = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt = fgmc('OU-NTS', [alphas(i), b, sigma, k, theta], N, M, T, 16, seed, 1e-8);
    [TCumulantsOUNTSalphaAsymm(i,:), ECumulantsOUNTSalphaAsymm(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'OU-NTS');
end

%% NTS-OU VARYING ALPHA, SYMMETRIC
alphas = [0.8, 0.6, 0.4, 0.2];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

TCumulantsNTSOUalphaSymm = zeros(numel(alphas), 4);
ECumulantsNTSOUalphaSymm = zeros(numel(alphas), 4);
timeNTSOUalphaSymm = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt = fgmc('NTS-OU', [alphas(i), b, sigma, k, theta], N, M, T, 24, seed, 1e-10);
    [TCumulantsNTSOUalphaSymm(i,:), ECumulantsNTSOUalphaSymm(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'NTS-OU');
end

%% NTS-OU VARYING ALPHA, ASYMMETRIC
alphas = [0.8, 0.6, 0.4, 0.2];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

TCumulantsNTSOUalphaAsymm = zeros(numel(alphas), 4);
ECumulantsNTSOUalphaAsymm = zeros(numel(alphas), 4);
timeNTSOUalphaAsymm = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt = fgmc('NTS-OU', [alphas(i), b, sigma, k, theta], N, M, T, 24, seed, 1e-10);
    [TCumulantsNTSOUalphaAsymm(i,:), ECumulantsNTSOUalphaAsymm(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'NTS-OU');
end

%% OU-TS VARYING ALPHA
alphas = [1.6, 1.2, 0.8, 0.4, -1.0, -2.0];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsOUTSalpha = zeros(numel(alphas), 4);
ECumulantsOUTSalpha = zeros(numel(alphas), 4);
timeOUTSalpha = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt= fgmc('OU-TS', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 18, seed, 1e-10);
    [TCumulantsOUTSalpha(i,:), ECumulantsOUTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS');
end

%% TS-OU VARYING ALPHA
alphas = [1.6, 1.2, 0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSOUalpha = zeros(numel(alphas), 4);
ECumulantsTSOUalpha = zeros(numel(alphas), 4);
timeTSOUalpha = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    Xt= fgmc('TS-OU', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 25, seed, 1e-12);
    [TCumulantsTSOUalpha(i,:), ECumulantsTSOUalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU');
end

%% Processes listed separately to run with a single alpha to better study the single case
% Note: the values of alpha here are to be intended as the last one used
% for a specific study during the developement of the code. 
% For the final results use the upper part of this script.

%% OU-NTS SYMMETRIC
alpha = -2; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

tic
Xt = fgmc('OU-NTS', [alpha, b, sigma, k, theta], N, M, T, 16, seed, 1e-8);
timeOUNTSsymm = toc;
[TCumulantsOUNTSsymm, ECumulantsOUNTSsymm] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'OU-NTS')

%% OU-NTS ASYMMETRIC
alpha = 0.4; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

tic
Xt = fgmc('OU-NTS', [alpha, b, sigma, k, theta], N, M, T, 16, seed, 0);
timeOUNTSasymm = toc;
[TCumulantsOUNTSasymm, ECumulantsOUNTSasymm] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'OU-NTS')

%% NTS-OU SYMMETRIC
alpha = 0.2; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0;

tic
Xt = fgmc('NTS-OU', [alpha, b, sigma, k, theta], N, M, T, 24, seed, 1e-10);
timeNTSOUsymm = toc;
[TCumulantsNTSOUsymm, ECumulantsNTSOUsymm] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'NTS-OU')

%% NTS-OU ASYMMETRIC
alpha = 0.2; b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

tic
Xt = fgmc('NTS-OU', [alpha, b, sigma, k, theta], N, M, T, 24, seed, 1e-10);
timeNTSOUasymm = toc;
[TCumulantsNTSOUasymm, ECumulantsNTSOUasymm] = computeCumulants(Xt(:,end), [alpha, b, sigma, k, theta], T, 'NTS-OU')

%% OU-TS
alpha = 0.4; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

tic
Xt = fgmc('OU-TS', [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 16, seed, 1e-10);
timeOUTS = toc;
[TCumulantsOUTS, ECumulantsOUTS] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS')

%% TS-OU
alpha = 0.8; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

tic
Xt = fgmc('TS-OU', [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, 25, seed, 1e-12);
timeTSOU = toc;
[TCumulantsTSOU, ECumulantsTSOU] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU')
