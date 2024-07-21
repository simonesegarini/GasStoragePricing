%% exDe - IMPLEMENTATION OF SABINO AND SABINO & CUFARO PETRONI 2022 PAPERS
clear all, close all, clc

%% SIMULATION PARAMETERS
T = 1/12; % Time horizon
M = 1; % Number of steps
N = 1e7; % Number of simulations
seed = 2; % Seed for the random sampling

%% OU-TS VARYING ALPHA
alphas = [0.8, 0.4, -1.0, -2.0];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsOUTSalpha = zeros(numel(alphas), 4);
ECumulantsOUTSalpha = zeros(numel(alphas), 4);
timeOUTSalpha = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    tic
    Xt = exDe('OU-TS', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed);
    timeOUTSalpha(i) = toc;
    [TCumulantsOUTSalpha(i,:), ECumulantsOUTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS');
end

%% TS-OU VARYING ALPHA
alphas = [0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSOUalpha = zeros(numel(alphas), 4);
ECumulantsTSOUalpha = zeros(numel(alphas), 4);
timeTSOUalpha = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    tic
    Xt = exDe('TS-OU', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed);
    timeTSOUalpha(i) = toc;
    [TCumulantsTSOUalpha(i,:), ECumulantsTSOUalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU');
end

%% OU-TS
alpha = 0.7; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

tic
Xt = exDe('OU-TS', [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed);
timeOUTS = toc;
[TCumulantsOUTS, ECumulantsOUTS] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS')

%% TS-OU
alpha = 0.7; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

tic
Xt = exDe('TS-OU', [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed);
timeTSOU = toc;
[TCumulantsTSOU, ECumulantsTSOU] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU')