clear all, close all, clc

T = 10;
S0 = 14.88;
seed = 2;

paramsTS = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case
% Number of simulations.
numberSimulations = 500; 

tic
X_OU_TS = spotSimulation('OU-TS', paramsTS, numberSimulations, T, T, 18, seed, 1e-10);
timeOU_TS = toc;
S_OU_TS = S0*exp(X_OU_TS);

% Simulation of the underlying for the TS-OU case.
tic
X_TS_OU = spotSimulation('TS-OU', paramsTS, numberSimulations, T, T, 25, seed, 1e-12, 'exde');
timeTS_OU = toc;
S_TS_OU = S0*exp(X_TS_OU);

simulations = [1, 10, 25, 50];

% OU-TS
plotSpot(simulations, S_OU_TS, 'OU-TS', 0)

% TS-OU
plotSpot(simulations, S_TS_OU, 'TS-OU', 0)

%%
[TCumulantsOUTS, ECumulantsOUTS] = computeCumulants(X_OU_TS(:,end), paramsTS, T, 'OU-TS')
[TCumulantsTSOU, ECumulantsTSOU] = computeCumulants(X_TS_OU(:,end), paramsTS, T, 'TS-OU')