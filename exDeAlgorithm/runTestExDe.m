%% exDe - IMPLEMENTATION OF SABINO AND SABINO & CUFARO PETRONI 2022 PAPERS
clear all, close all, clc

%% SIMULATION PARAMETERS
T = 1; % Time horizon
M = 1; % Number of steps
N = 1e7; % Number of simulations
seed = 4; % Seed for the random sampling
%% OU-TS
alpha = 0.7; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

tic
Xt = exDe('OU-TS', [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], N, M, T, seed, 'SSR');
timeOUTS = toc;
[TCumulantsOUTS, ECumulantsOUTS] = computeCumulants(Xt(:,end), [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], T, 'OU-TS')

%%

% figure
% plot([0, T], Xt(1:10, :), 'LineWidth', 3)