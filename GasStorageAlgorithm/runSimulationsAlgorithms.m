clear all, close all, clc

%% Set up FLAGS for processes.
OU_TS_FGMC = 21; TS_OU_FGMC = 31; 
OU_TS_SSR = 22; TS_OU_SSR = 32;
OU_TS_DR = 23; TS_OU_DR = 33;

OU_NTS = 4; NTS_OU = 5;

%% Define the struct for processes simulations.

simulation.Simulations = 1e7; % Number of paths simulated.
simulation.Maturity = 1/12; % Maturity expressed in years.
simulation.Steps = 1; % Steps for simulation.
simulation.Seed = 2; % Seed for the rng, NOT mandatory.

% SPOT PRICES SIMULATIONS.
% Process available: OU, OU-NTS, NTS-OU, OU-TS, TS-OU with also all special
% cases.
% If OU-NTS or NTS-OU are selected, parameters should be given as 
% [alpha, b, sigma, k, theta].
% If OU-TS or TS-OU are selected, parameters should be given as 
% [b, beta_p, beta_n, c_p, c_n, gamma_c].
%
% Stability parameters alphas are to be given in another vector so the
% script can run for all the desired values.
simulation.Process = OU_TS_SSR;
simulation.Parameters = [0.1, 2.5, 3.5, 0.5, 1, 0]; % Parameters for TS.
% simulation.Parameters = [0.2162, 0.201, 0.256, 0]; % Parameters for NTS.
simulation.Alphas = [0.8, 0.4, -1.0, -2.0];
% simulation.Alphas = [0.8];

%% Run the algorithm

[simulation.X, simulation.ECumulants, simulation.TCumulants, simulation.times] ...
    = simulationAlgorithms(simulation);

displaySimulationResults(simulation)
%% Save the environment
save('simulation_environment.mat')