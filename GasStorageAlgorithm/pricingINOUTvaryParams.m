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

% Setted for OU-NTS, both sym and asym.
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0]; thetas = [0, 0.1];

% Number of simulations.
numberSimulations = [250, 2500, 25000];

% Discretization for the backward induction.
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

% Choose the model.
model = 'OU-NTS';

% Allocate the values to be saved.
price_IN = zeros(3,1); IN_STD = zeros(3,1); price_IN_CI = zeros(3,2);
price_OUT = zeros(3,1); OUT_STD = zeros(3,1); price_OUT_CI = zeros(3,2);

for i=1:numel(thetas)
    for j=1:numel(alphas)
        for z=1:numel(numberSimulations)

            numSims = numberSimulations(z);

            % Params for the OU-NTS case.
            params = [alphas(j), 0.2162, 0.201, 0.256, thetas(i)]; 

            % Simulation of the underlying for the OU-NTS case
            X = spotSimulation(model, params, numSims, 365, T, 16, 1, 1e-8);
            S = S0*exp(X);
    
            % 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
            % Matrices with the cashflows
            cashflows_IN = penFunc(S(:,end), ones(numSims,1)*dV');
            
            % 3. APPLY BACKWARD INDUCTION
            [cashflows_IN, ~, regressors_IN] = backwardInduction(S, cashflows_IN, payoff, N, numSims, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

            X_OUT = spotSimulation(model, params, numSims, 365, T, 16, 2, 1e-8);
            S_OUT = S0*exp(X_OUT);
            cashflows_OUT = penFunc(S_OUT(:,end), ones(numSims,1)*dV');
            cashflows_OUT = priceOut(S_OUT, cashflows_OUT, regressors_IN, payoff, N, numSims, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
            
            % 4. PRICE
            % Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
            [price_IN(z,:), IN_STD(z,:), price_IN_CI(z,:)] = normfit(cashflows_IN(:, index_V0));
            [price_OUT(z,:), OUT_STD(z,:), price_OUT_CI(z,:)] = normfit(cashflows_OUT(:, index_V0));
        end
        name = sprintf('INOUT_Model_%s__alpha_%0.1f__theta_%0.1f.mat', model, alphas(j), thetas(i));
        save(name, 'price_IN', 'IN_STD', 'price_IN_CI', 'price_OUT', 'OUT_STD', 'price_OUT_CI')
    end
end