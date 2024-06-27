clear all, close all, clc

% This script is used to simulate prices evolution and to price the
% contract by varying parameters when more are proposed in papers.
% Here is priced with OU process. Results are saved in .mat files.

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

% volume constraints and values.
Vmin = 0;
Vmax = 250000;
V0 = 100000;
VT = 100000;

% deltaV variations constraints
Imin = @(v) max(Vmin - v, -7500);
Imax = @(v) min(Vmax - v, 2500);

% penalization functions
penFunc = @(s, v) -s.*abs(v-VT).^2;

%% 1. SPOT PRICE SIMULATION 
T = 365;
S0 = 14.88;

% OU-NTS
vols = [0.0315, 0.0945];

% Number of simulations
numberSimulations = 500; 

% Discretization for the backward induction
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

model = 'OU';

for j=1:numel(vols)
    paramsOU = [vols(j), 0.05, 0]; % Parameters for the OU case
    
    % Simulation of the underlying for the OU-NTS case
    Xs = spotSimulation(model, paramsOU, numberSimulations, 365, T, 0, 1, 0);
    X_OU = Xs(1:numberSimulations, :); XAV_OU = Xs(numberSimulations+1:end, :);
    S_OU = S0*exp(X_OU);
    SAV_OU = S0*exp(XAV_OU);


    % 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
    % Matrices with the cashflows
    cashflows_OU_T = penFunc(S_OU(:,end), ones(numberSimulations,1)*dV');
    cashflows_OU_AV_T = penFunc(SAV_OU(:,end), ones(numberSimulations,1)*dV');
    
    % 3. APPLY BACKWARD INDUCTION
    cashflows_OU = priceIn(S_OU, cashflows_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
    cashflows_OU_AV = priceIn(SAV_OU, cashflows_OU_AV_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

    % 4. PRICE
    % Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
    [price_OU_IN, IN_OU_STD, price_OU_IN_CI] = normfit(0.5*(cashflows_OU(:,index_V0) + cashflows_OU_AV(:,index_V0)));

    name = sprintf('Model_%s__vol_%0.4f.mat', model, vols(j));
    save(name, 'S_OU', 'SAV_OU', 'price_OU_IN', 'IN_OU_STD', 'price_OU_IN_CI')
end
