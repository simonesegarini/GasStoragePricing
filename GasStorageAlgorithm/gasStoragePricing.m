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

% OU.
paramsOUL = [0.0315, 0.05, 0]; % Parameters for the OU case, low volatility.
paramsOUH = [0.0945, 0.05, 0]; % Parameters for the OU case, high volatility.

% OU-NTS.
paramsNTSsymm = [0.6, 0.2162, 0.201, 0.256, 0]; % Params for the OU-NTS symmetric case.
paramsNTSasymm = [0.6, 0.2162, 0.201, 0.256, 0.1]; % Params for the OU-NTS asymmetric case.

% OU-TS.
paramsTS = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case

% Number of simulations.
numberSimulations = 500; 

% Discretization for the backward induction.
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

%% Simulations with all type of processes studied in fgmc algorithm
% Simulation of the underlying for the OU case, low volatility.
XsL = spotSimulation('OU', paramsOUL, numberSimulations, 365, T, 0, seed, 0);
X_OUL = XsL(1:numberSimulations, :); XAV_OUL = XsL(numberSimulations+1:end, :);
S_OUL = S0*exp(X_OUL);
SAV_OUL = S0*exp(XAV_OUL);

% Simulation of the underlying for the OU case, high volatility.
XsH = spotSimulation('OU', paramsOUH, numberSimulations, 365, T, 0, seed, 0);
X_OUH = XsH(1:numberSimulations, :); XAV_OUH = XsH(numberSimulations+1:end, :);
S_OUH = S0*exp(X_OUH);
SAV_OUH = S0*exp(XAV_OUH);

% Simulation of the underlying for the OU-NTS symmetric case.
tic
X_OU_NTSsymm = spotSimulation('OU-NTS', paramsNTSsymm, numberSimulations, 365, T, 16, seed, 1e-8);
timeOU_NTSsymm = toc;
S_OU_NTSsymm = S0*exp(X_OU_NTSsymm);

% Simulation of the underlying for the OU-NTS asymmetric case.
tic
X_OU_NTSasymm = spotSimulation('OU-NTS', paramsNTSasymm, numberSimulations, 365, T, 16, seed, 1e-8);
timeOU_NTSasymm = toc;
S_OU_NTSasymm = S0*exp(X_OU_NTSasymm);

% Simulation of the underlying for the NTS-OU symmetric case.
tic
X_NTS_OUsymm = spotSimulation('NTS-OU', paramsNTSsymm, numberSimulations, 365, T, 24, seed, 1e-10);
timeNTS_OUsymm = toc;
S_NTS_OUsymm = S0*exp(X_NTS_OUsymm);

% Simulation of the underlying for the NTS-OU asymmetric case.
tic
X_NTS_OUasymm = spotSimulation('NTS-OU', paramsNTSasymm, numberSimulations, 365, T, 24, seed, 1e-10);
timeNTS_OUasymm = toc;
S_NTS_OUasymm = S0*exp(X_NTS_OUasymm);

% Simulation of the underlying for the OU-TS case.
tic
X_OU_TS = spotSimulation('OU-TS', paramsTS, numberSimulations, 365, T, 18, seed, 1e-10, 'exde');
timeOU_TS = toc;
S_OU_TS = S0*exp(X_OU_TS);

% Simulation of the underlying for the TS-OU case.
tic
X_TS_OU = spotSimulation('TS-OU', paramsTS, numberSimulations, 365, T, 25, seed, 1e-12, 'exde');
timeTS_OU = toc;
S_TS_OU = S0*exp(X_TS_OU);

%% Plot to see some paths
simulations = [1, 5, 25, 50];

% OU, low and high volatilities, normal and AV
plotSpot(simulations, S_OUL, 'OU-L', 1)
plotSpot(simulations, S_OUH, 'OU-H', 1)
plotSpot(simulations, SAV_OUL, 'OU-AVL', 1)
plotSpot(simulations, SAV_OUH, 'OU-AVH', 1)

% OU-NTS, symmetric and asymmetric
plotSpot(simulations, S_OU_NTSsymm, 'OU-NTSsymm', 1)
plotSpot(simulations, S_NTS_OUsymm, 'NTS-OUsymm', 1)

% NTS-OU symmetric and asymmetric
plotSpot(simulations, S_OU_NTSasymm, 'OU-NTSasymm', 1)
plotSpot(simulations, S_NTS_OUasymm, 'NTS-OUasymm', 1)

% OU-TS
plotSpot(simulations, S_OU_TS, 'OU-TS', 1)

% TS-OU
plotSpot(simulations, S_TS_OU, 'TS-OU', 1)

%% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
% Assign the final penalization cost to the cashflows matrices, they have a
% T in the end to refer as the final time step, created this variable to
% avoid new initialization each time we want to try a different approach.

% OU, low and high volatilities, normal and AV
cashflows_OUL_T = penFunc(S_OUL(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVL_T = penFunc(SAV_OUL(:,end), ones(numberSimulations,1)*dV');
cashflows_OUH_T = penFunc(S_OUH(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVH_T = penFunc(SAV_OUH(:,end), ones(numberSimulations,1)*dV');

% OU-NTS, symmetric and asymmetric
cashflows_OU_NTSsymm_T = penFunc(S_OU_NTSsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_NTSasymm_T = penFunc(S_OU_NTSasymm(:,end), ones(numberSimulations,1)*dV');

% NTS-OU symmetric and asymmetric
cashflows_NTS_OUsymm_T = penFunc(S_NTS_OUsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUasymm_T = penFunc(S_NTS_OUasymm(:,end), ones(numberSimulations,1)*dV');

% OU-TS
cashflows_OU_TS_T = penFunc(S_OU_TS(:,end), ones(numberSimulations,1)*dV');

% TS-OU
cashflows_TS_OU_T = penFunc(S_TS_OU(:,end), ones(numberSimulations,1)*dV');

%% 3. APPLY BACKWARD INDUCTION WITH PRICEIN
% Use IN-pricing to avoid using too much space for the matrices of policies
% and for the regressors, they are needed just for plot and OUT-price
% purpose.

% OU, low and high volatilities, normal and AV
cashflows_OUL = priceIn(S_OUL, cashflows_OUL_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVL = priceIn(SAV_OUL, cashflows_OU_AVL_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OUH = priceIn(S_OUH, cashflows_OUH_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVH = priceIn(SAV_OUH, cashflows_OU_AVH_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-NTS, symmetric and asymmetric
cashflows_OU_NTSsymm = priceIn(S_OU_NTSsymm, cashflows_OU_NTSsymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_NTSasymm = priceIn(S_OU_NTSasymm, cashflows_OU_NTSasymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% NTS-OU symmetric and asymmetric
cashflows_NTS_OUsymm = priceIn(S_NTS_OUsymm, cashflows_NTS_OUsymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OUasymm = priceIn(S_NTS_OUasymm, cashflows_NTS_OUasymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-TS
cashflows_OU_TS = priceIn(S_OU_TS, cashflows_OU_TS_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

%TS-OU
cashflows_TS_OU = priceIn(S_TS_OU, cashflows_TS_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

%% 4. PRICE
% Compute the final price as the mean of accumulated cash flows at t=0
% across all simulations. We take prices and CIs from normfit while the
% std error from the function created for this.

% OU low volatility
[price_OUL_IN, ~, price_OUL_IN_CI] = normfit(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)));
IN_OUL_STD = standardError(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)));

% OU high volatility
[price_OUH_IN, ~, price_OUH_IN_CI] = normfit(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)));
IN_OUH_STD = standardError(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)));

% OU-NTS symmetric
[price_OU_NTSsymm_IN, ~, price_OU_NTSsymm_IN_CI] = normfit(cashflows_OU_NTSsymm(:, index_V0));
IN_OU_NTSsymm_STD = standardError(cashflows_OU_NTSsymm(:, index_V0));

% OU-NTS asymmetric
[price_OU_NTSasymm_IN, ~, price_OU_NTSasymm_IN_CI] = normfit(cashflows_OU_NTSasymm(:, index_V0));
IN_OU_NTSasymm_STD = standardError(cashflows_OU_NTSasymm(:, index_V0));

% NTS-OU symmetric
[price_NTS_OUsymm_IN, ~, price_NTS_OUsymm_IN_CI] = normfit(cashflows_NTS_OUsymm(:, index_V0));
IN_NTS_OUsymm_STD = standardError(cashflows_NTS_OUsymm(:, index_V0));

% NTS-OU asymmetric
[price_NTS_OUasymm_IN, ~, price_NTS_OUasymm_IN_CI] = normfit(cashflows_NTS_OUasymm(:, index_V0));
IN_NTS_OUasymm_STD = standardError(cashflows_NTS_OUasymm(:, index_V0));

% OU-TS
[price_OU_TS_IN, ~, price_OU_TS_IN_CI] = normfit(cashflows_OU_TS(:, index_V0));
IN_OU_TS_STD = standardError(cashflows_OU_TS(:, index_V0));

% TS-OU
[price_TS_OU_IN, ~, price_TS_OU_IN_CI] = normfit(cashflows_TS_OU(:, index_V0));
IN_TS_OU_STD = standardError(cashflows_TS_OU(:, index_V0));

%% Plots of distribution of storage prices
% Plot the distribution of the IN price (save the plots if the last inputs is 1)

% OU case, low volatility
plotPriceDistribution(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)), numberSimulations, 'OU-L', paramsOUL(1)*100, 1);

% OU case, high volatility
plotPriceDistribution(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)), numberSimulations, 'OU-H', paramsOUH(1)*100, 1);

% OU-NTS symmetric case
plotPriceDistribution(cashflows_OU_NTSsymm(:, index_V0), numberSimulations, 'OU-NTSsymm', paramsNTSsymm(1), 1);

% OU-NTS asymmetric case
plotPriceDistribution(cashflows_OU_NTSasymm(:, index_V0), numberSimulations, 'OU-NTSasymm', paramsNTSasymm(1), 1);

% NTS-OU symmetric case
plotPriceDistribution(cashflows_NTS_OUsymm(:, index_V0), numberSimulations, 'NTS-OUsymm', paramsNTSsymm(1), 1);

% NTS-OU asymmetric case
plotPriceDistribution(cashflows_NTS_OUasymm(:, index_V0), numberSimulations, 'NTS-OUasymm', paramsNTSasymm(1), 1);

% OU-TS case
plotPriceDistribution(cashflows_OU_TS(:, index_V0), numberSimulations, 'OU-TS', paramsTS(1), 1);

% TS-OU case
plotPriceDistribution(cashflows_TS_OU(:, index_V0), numberSimulations, 'TS-OU', paramsTS(1), 1);