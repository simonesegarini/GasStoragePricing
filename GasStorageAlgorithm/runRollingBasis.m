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
paramsNTSsymm = [0.6, 0.2162, 0.201, 0.256, 0]; % Params for the OU-NTS symmetric case.
paramsNTSasymm = [0.6, 0.2162, 0.201, 0.256, 0.1]; % Params for the OU-NTS asymmetric case.

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
meansSym = zeros(numel(order), numel(numberSimulations));
stdErrorsSym = zeros(numel(order), numel(numberSimulations));
timesSym = zeros(numel(order), numel(numberSimulations));
CIsSym = zeros(numel(order), numel(numberSimulations), 2);

meansAsym = zeros(numel(order), numel(numberSimulations));
stdErrorsAsym = zeros(numel(order), numel(numberSimulations));
timesAsym = zeros(numel(order), numel(numberSimulations));
CIsAsym = zeros(numel(order), numel(numberSimulations), 2);


for simIt=1:numel(numberSimulations)
    
    
    Xsym = spotSimulation('NTS-OU', paramsNTSsymm, numberSimulations(simIt), 365, T, 24, seed, 1e-10);
    Ssym = S0*exp(Xsym);
    Xasym = spotSimulation('NTS-OU', paramsNTSasymm, numberSimulations(simIt), 365, T, 24, seed, 1e-10);
    Sasym = S0*exp(Xasym);

    for regIt=1:numel(order)
        cashflowsSym = penFunc(Ssym(:,end), ones(numberSimulations(simIt),1)*dV');
        cashflowsAsym = penFunc(Sasym(:,end), ones(numberSimulations(simIt),1)*dV');
        
        tic
        cashflowsSym = priceIn(Ssym, cashflowsSym, payoff, N, numberSimulations(simIt), delta, alpha, T, maxInjection, maxWithdraw, order(regIt), type(regIt));
        timesSym(regIt, simIt) = toc;

        tic
        cashflowsAsym = priceIn(Sasym, cashflowsAsym, payoff, N, numberSimulations(simIt), delta, alpha, T, maxInjection, maxWithdraw, order(regIt), type(regIt));
        timesAsym(regIt, simIt) = toc;
    
        [meansSym(regIt, simIt), ~, CIsSym(regIt, simIt, :)] = normfit(cashflowsSym(:,index_V0));
        stdErrorsSym(regIt, simIt) = standardError(cashflowsSym(:,index_V0));

        [meansAsym(regIt, simIt), ~, CIsAsym(regIt, simIt, :)] = normfit(cashflowsAsym(:,index_V0));
        stdErrorsAsym(regIt, simIt) = standardError(cashflowsAsym(:,index_V0));

    end
end
