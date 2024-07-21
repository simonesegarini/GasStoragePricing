clear all, close all, clc

%% STANDARD STORAGE CONTRACT DATA
% Cost of injection parameters.
a1 = 0; b1 = 0;

% Cost of selling parameters.
a2 = 0; b2 = 0;

% Discount rate.
delta = 0;

% Payoff for each action.
payoff = @(s, deltaV) -((1+a1).*s + b1).*deltaV.*(deltaV > 0) ...
                - ((1-a2).*s - b2).*deltaV.*(deltaV < 0);

% Volume constraints and values.
Vmin = 0;
Vmax = 250000;
V0 = 100000;
VT = 100000;

% Maximum and minimum allowed injection and withdrawal of gas in time.
Imin = @(v) max(Vmin - v, -7500);
Imax = @(v) min(Vmax - v, 2500);

% Definition of penalization function in order to force the final desired
% storage volume.
penFunc = @(s, v) -s.*abs(v-VT).^2;

%% 1. SPOT PRICE SIMULATION 
T = 365;
S0 = 14.88;
seed = 2;

% Parameters.
% OU-TS.
paramsTS = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case

% Number of simulations.
numberSimulations = [50, 500, 5000]; 

% Discretization for the backward induction.
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

% Define vector for the basis we want to try
type = ["polynomial", "polynomial", "polynomial", "bspline", "bspline"];
order = [3,4,5,5,6];

% Define variables to store data needed for comparison.
means = zeros(numel(order), numel(numberSimulations));
stdErrors = zeros(numel(order), numel(numberSimulations));
times = zeros(numel(order), numel(numberSimulations));
CIs = zeros(numel(order), numel(numberSimulations), 2);


for simIt=1:numel(numberSimulations)
    
    
    X = spotSimulation('OU-TS', paramsTS, numberSimulations(simIt), 365, T, 24, seed, 1e-10, 'exde');
    S = S0*exp(X);

    for regIt=1:numel(order)
        cashflowsSym = penFunc(S(:,end), ones(numberSimulations(simIt),1)*dV');
        
        tic
        cashflowsSym = priceIn(S, cashflowsSym, payoff, N, numberSimulations(simIt), delta, alpha, T, maxInjection, maxWithdraw, order(regIt), type(regIt));
        times(regIt, simIt) = toc;
    
        [means(regIt, simIt), ~, CIs(regIt, simIt, :)] = normfit(cashflowsSym(:,index_V0));
        stdErrors(regIt, simIt) = standardError(cashflowsSym(:,index_V0));


    end
end
