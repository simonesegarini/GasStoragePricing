clear all, close all, clc

% TO ADJUST WITH MODIFICATIONS

%% STANDARD STORAGE CONTRACT DATA
% cost of injection
a1 = 0; b1 = 0;

% cost of selling
a2 = 0; b2 = 0;

% discount rate
delta = 0;

% payoff
h = @(s, deltaV) -((1+a1).*s + b1).*deltaV.*(deltaV > 0) ...
                - ((1-a2).*s - b2).*deltaV.*(deltaV < 0);

% time reference
ttm = 1;
startDate= datetime(2005, 07, 01);
endDate = datetime(2006, 06, 30);

% volume constraints and values
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

% OU-TS
alphas = [1.6, 1.2, 0.8, 0.4];

% Number of simulations
M = 500; 

% Discretization for the backward induction
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

model = 'OU-TS';

for j=1:numel(alphas)
    paramsTS = [alphas(j), 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case
    
    % Simulation of the underlying for the OU-NTS case
    X_TS = spotSimulation(model, paramsTS, M, 365, T, 16, 1);
    S_TS = S0*exp(X_TS);


    % 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
    % Matrices with the cashflows
    cashflows_TS = penFunc(S_TS(:,end), ones(M,1)*dV');
    
    % 3. APPLY BACKWARD INDUCTION
    cashflows_TS = priceIn(S_TS, cashflows_TS, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);

    % 4. PRICE
    % Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
    [price_TS_IN, IN_TS_STD, price_TS_IN_CI] = normfit(cashflows_TS(:, index_V0));

    name = sprintf('Model_%s__alpha_%0.1f.mat', model, alphas(j));
    save(name, 'S_TS', 'price_TS_IN', 'IN_TS_STD', 'price_TS_IN_CI')
end
