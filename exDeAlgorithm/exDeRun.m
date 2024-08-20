%% exDe - IMPLEMENTATION OF SABINO AND SABINO & CUFARO PETRONI 2022 PAPERS
clear all, close all, clc

%% SIMULATION PARAMETERS
T = 1; % Time horizon
M = 1; % Number of steps
N = 1e7; % Number of simulations
seed = 2; % Seed for the random sampling

%% OU-TS VARYING ALPHA SSR
% alphas = [0.8, 0.4, -1.0, -2.0];
% alphas = [0.4, -1.0, -2.0];
alphas = [0.8];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsOUTSalphaSSR = zeros(numel(alphas), 4);
ECumulantsOUTSalphaSSR = zeros(numel(alphas), 4);
timeOUTSalphaSSR = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    disp(['OU-TS SSR with \alpha = ', num2str(alphas(i))])
    tic
    Xt = exDe('OU-TS', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed, 'SSR');
    timeOUTSalphaSSR(i) = toc;
    disp(['Time for OU-TS SSR with \alpha = ', num2str(alphas(i)), ' is ', num2str(timeOUTSalphaSSR(i))])
    [TCumulantsOUTSalphaSSR(i,:), ECumulantsOUTSalphaSSR(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS');
end

%% OU-TS VARYING ALPHA DR
alphas = [0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsOUTSalphaDR = zeros(numel(alphas), 4);
ECumulantsOUTSalphaDR = zeros(numel(alphas), 4);
timeOUTSalphaDR = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    disp(['OU-TS DR with \alpha = ', num2str(alphas(i))])
    tic
    Xt = exDe('OU-TS', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed, 'DR');
    timeOUTSalphaDR(i) = toc;
    disp(['Time for OU-TS DR with \alpha = ', num2str(alphas(i)), ' is ', num2str(timeOUTSalphaDR(i))])
    [TCumulantsOUTSalphaDR(i,:), ECumulantsOUTSalphaDR(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS');
end

%% TS-OU VARYING ALPHA SSR
alphas = [0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSOUalphaSSR = zeros(numel(alphas), 4);
ECumulantsTSOUalphaSSR = zeros(numel(alphas), 4);
timeTSOUalphaSSR = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    disp(['TS-OU SSR with \alpha = ', num2str(alphas(i))])
    tic
    Xt = exDe('TS-OU', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed, 'SSR');
    timeTSOUalphaSSR(i) = toc;
    disp(['Time for TS-OU SSR with \alpha = ', num2str(alphas(i)), ' is ', num2str(timeTSOUalphaSSR(i))])
    [TCumulantsTSOUalphaSSR(i,:), ECumulantsTSOUalphaSSR(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU');
end

%% TS-OU VARYING ALPHA DR
alphas = [0.8, 0.4];
b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

TCumulantsTSOUalphaDR = zeros(numel(alphas), 4);
ECumulantsTSOUalphaDR = zeros(numel(alphas), 4);
timeTSOUalphaDR = zeros(numel(alphas), 1);

for i=1:numel(alphas)
    disp(['TS-OU DR with \alpha = ', num2str(alphas(i))])
    tic
    Xt = exDe('TS-OU', [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed, 'DR');
    timeTSOUalphaDR(i) = toc;
    disp(['Time for TS-OU DR with \alpha = ', num2str(alphas(i)), ' is ', num2str(timeTSOUalphaDR(i))])
    [TCumulantsTSOUalphaDR(i,:), ECumulantsTSOUalphaDR(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'TS-OU');
end