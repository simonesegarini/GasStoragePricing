clear all, close all, clc

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

% OU-NTS
alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0]; thetas = [0, 0.1];

% Number of simulations
Ms = [250, 2500, 25000];

% Discretization for the backward induction
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

model = 'OU-NTS';

price_NTS_IN = zeros(3,1); IN_NTS_STD = zeros(3,1); price_NTS_IN_CI = zeros(3,2);
price_NTS_OUT = zeros(3,1); OUT_NTS_STD = zeros(3,1); price_NTS_OUT_CI = zeros(3,2);

for i=1:numel(thetas)
    for j=1:numel(alphas)
        for z=1:numel(Ms)

            M = Ms(z);

            % OU-NTS
            paramsNTS = [alphas(j), 0.2162, 0.201, 0.256, thetas(i)]; % Params for the OU-NTS case

            % Simulation of the underlying for the OU-NTS case
            X_NTS = spotSimulation(model, paramsNTS, M, 365, T, 16, 1);
            S_NTS = S0*exp(X_NTS);
    
            % 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
            % Matrices with the cashflows
            cashflows_NTS = penFunc(S_NTS(:,end), ones(M,1)*dV');
            
            % 3. APPLY BACKWARD INDUCTION
            [cashflows_NTS, policies_NTS, regression_NTS] = backwardInduction(S_NTS, cashflows_NTS, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);

            X_NTS_OUT = spotSimulation(model, paramsNTS, M, 365, T, 16, 2);
            S_NTS_OUT = S0*exp(X_NTS_OUT);
            cashflows_NTS_OUT = penFunc(S_NTS_OUT(:,end), ones(M,1)*dV');
            cashflows_NTS_OUT = priceOut(S_NTS_OUT, cashflows_NTS_OUT, regression_NTS, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
            
            % 4. PRICE
            % Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
            [price_NTS_IN(z,:), IN_NTS_STD(z,:), price_NTS_IN_CI(z,:)] = normfit(cashflows_NTS(:, index_V0));
            [price_NTS_OUT(z,:), OUT_NTS_STD(z,:), price_NTS_OUT_CI(z,:)] = normfit(cashflows_NTS_OUT(:, index_V0));
        end
        name = sprintf('INOUT_Model_%s__alpha_%0.1f__theta_%0.1f.mat', model, alphas(j), thetas(i));
        save(name, 'price_NTS_IN', 'IN_NTS_STD', 'price_NTS_IN_CI', 'price_NTS_OUT', 'OUT_NTS_STD', 'price_NTS_OUT_CI')
    end
end