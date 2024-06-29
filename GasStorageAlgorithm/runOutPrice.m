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
X_OU_NTSsymm = spotSimulation('OU-NTS', paramsNTSsymm, numberSimulations, 365, T, 16, seed, 1e-8);
S_OU_NTSsymm = S0*exp(X_OU_NTSsymm);

% Simulation of the underlying for the OU-NTS asymmetric case.
X_OU_NTSasymm = spotSimulation('OU-NTS', paramsNTSasymm, numberSimulations, 365, T, 16, seed, 1e-8);
S_OU_NTSasymm = S0*exp(X_OU_NTSasymm);

% Simulation of the underlying for the NTS-OU symmetric case.
X_NTS_OUsymm = spotSimulation('NTS-OU', paramsNTSsymm, numberSimulations, 365, T, 24, seed, 1e-10);
S_NTS_OUsymm = S0*exp(X_NTS_OUsymm);

% Simulation of the underlying for the NTS-OU asymmetric case.
X_NTS_OUasymm = spotSimulation('NTS-OU', paramsNTSasymm, numberSimulations, 365, T, 24, seed, 1e-10);
S_NTS_OUasymm = S0*exp(X_NTS_OUasymm);

% Simulation of the underlying for the OU-TS case.
X_OU_TS = spotSimulation('OU-TS', paramsTS, numberSimulations, 365, T, 18, seed, 1e-10);
S_OU_TS = S0*exp(X_OU_TS);

% Simulation of the underlying for the TS-OU case.
X_TS_OU = spotSimulation('TS-OU', paramsTS, numberSimulations, 365, T, 25, seed, 1e-12);
S_TS_OU = S0*exp(X_TS_OU);

%% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
% Assign the final penalization cost to the cashflows matrices, they have a
% T in the end to refer as the final time step, created this variable to
% avoid new initialization each time we want to try a different approach.

% OU, low and high volatilities, normal and AV.
cashflows_OUL_T = penFunc(S_OUL(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVL_T = penFunc(SAV_OUL(:,end), ones(numberSimulations,1)*dV');
cashflows_OUH_T = penFunc(S_OUH(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVH_T = penFunc(SAV_OUH(:,end), ones(numberSimulations,1)*dV');

% OU-NTS, symmetric and asymmetric.
cashflows_OU_NTSsymm_T = penFunc(S_OU_NTSsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_NTSasymm_T = penFunc(S_OU_NTSasymm(:,end), ones(numberSimulations,1)*dV');

% NTS-OU symmetric and asymmetric.
cashflows_NTS_OUsymm_T = penFunc(S_NTS_OUsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUasymm_T = penFunc(S_NTS_OUasymm(:,end), ones(numberSimulations,1)*dV');

% OU-TS.
cashflows_OU_TS_T = penFunc(S_OU_TS(:,end), ones(numberSimulations,1)*dV');

% TS-OU.
cashflows_TS_OU_T = penFunc(S_TS_OU(:,end), ones(numberSimulations,1)*dV');

%% 3. APPLY BACKWARD INDUCTION
% OUT-price, run the backwardInduction algorithm that gives
% the regressors and the policies, run the part of priceOut after
% simulating the prices with another seed. (For now available only for 
% polynomial regression)
% After the computation of the regressors and the policies we can plot the
% volume evolution through time and the price using the out method.

% Get regressors, priceIn and policies.
% OU, low volatility
[cashflows_OUL, policies_OUL, regressors_OUL] = backwardInduction(S_OUL, cashflows_OUL_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
[cashflows_OU_AVL, policies_OU_AVL, regressors_OU_AVL] = backwardInduction(SAV_OUL, cashflows_OU_AVL_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU, high volatility.
[cashflows_OUH, policies_OUH, regressors_OUH] = backwardInduction(S_OUH, cashflows_OUH_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
[cashflows_OU_AVH, policies_OU_AVH, regressors_OU_AVH] = backwardInduction(SAV_OUH, cashflows_OU_AVH_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-NTS symmetric.
[cashflows_OU_NTSsymm, policies_OU_NTSsymm, regressors_OU_NTSsymm] = backwardInduction(S_OU_NTSsymm, cashflows_OU_NTSsymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-NTS asymmetric.
[cashflows_OU_NTSasymm, policies_OU_NTSasymm, regressors_OU_NTSasymm] = backwardInduction(S_OU_NTSasymm, cashflows_OU_NTSasymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% NTS-OU symmetric.
[cashflows_NTS_OUsymm, policies_NTS_OUsymm, regressors_NTS_OUsymm] = backwardInduction(S_NTS_OUsymm, cashflows_NTS_OUsymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% NTS-OU asymmetric.
[cashflows_NTS_OUasymm, policies_NTS_OUasymm, regressors_NTS_OUasymm] = backwardInduction(S_NTS_OUasymm, cashflows_NTS_OUasymm_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-TS.
[cashflows_OU_TS, policies_OU_TS, regressors_OU_TS] = backwardInduction(S_OU_TS, cashflows_OU_TS_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% TS-OU.
[cashflows_TS_OU, policies_TS_OU, regressors_TS_OU] = backwardInduction(S_TS_OU, cashflows_TS_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

%% 4. PRICEIN, TO HAVE A STARTING VALUE TO MATCH WITH PRICEOUT
% OU low volatility.
[price_OUL_IN, ~, price_OUL_IN_CI] = normfit(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)));
IN_OUL_STD = standardError(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)));

% OU high volatility.
[price_OUH_IN, ~, price_OUH_IN_CI] = normfit(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)));
IN_OUH_STD = standardError(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)));

% OU-NTS symmetric.
[price_OU_NTSsymm_IN, ~, price_OU_NTSsymm_IN_CI] = normfit(cashflows_OU_NTSsymm(:, index_V0));
IN_OU_NTSsymm_STD = standardError(cashflows_OU_NTSsymm(:, index_V0));

% OU-NTS asymmetric.
[price_OU_NTSasymm_IN, ~, price_OU_NTSasymm_IN_CI] = normfit(cashflows_OU_NTSasymm(:, index_V0));
IN_OU_NTSasymm_STD = standardError(cashflows_OU_NTSasymm(:, index_V0));

% NTS-OU symmetric.
[price_NTS_OUsymm_IN, ~, price_NTS_OUsymm_IN_CI] = normfit(cashflows_NTS_OUsymm(:, index_V0));
IN_NTS_OUsymm_STD = standardError(cashflows_NTS_OUsymm(:, index_V0));

% NTS-OU asymmetric.
[price_NTS_OUasymm_IN, ~, price_NTS_OUasymm_IN_CI] = normfit(cashflows_NTS_OUasymm(:, index_V0));
IN_NTS_OUasymm_STD = standardError(cashflows_NTS_OUasymm(:, index_V0));

% OU-TS.
[price_OU_TS_IN, ~, price_OU_TS_IN_CI] = normfit(cashflows_OU_TS(:, index_V0));
IN_OU_TS_STD = standardError(cashflows_OU_TS(:, index_V0));

% TS-OU.
[price_TS_OU_IN, ~, price_TS_OU_IN_CI] = normfit(cashflows_TS_OU(:, index_V0));
IN_TS_OU_STD = standardError(cashflows_TS_OU(:, index_V0));

%% Plot the evolution of the volume.
volumePlot(row_indices, T, V0, policies_OUL, alpha, 'OU-L', index_V0);
volumePlot(row_indices, T, V0, policies_OU_AVL, alpha, 'OU-L AV', index_V0);
volumePlot(row_indices, T, V0, policies_OUH, alpha, 'OU-H', index_V0);
volumePlot(row_indices, T, V0, policies_OU_AVH, alpha, 'OU-H AV', index_V0);
volumePlot(row_indices, T, V0, policies_OU_NTSsymm, alpha, 'OU-NTS sym', index_V0);
volumePlot(row_indices, T, V0, policies_OU_NTSasymm, alpha, 'OU-NTS asym', index_V0);
volumePlot(row_indices, T, V0, policies_NTS_OUsymm, alpha, 'NTS-OU sym', index_V0);
volumePlot(row_indices, T, V0, policies_NTS_OUasymm, alpha, 'NTS-OU asym', index_V0);
volumePlot(row_indices, T, V0, policies_OU_TS, alpha, 'OU-TS', index_V0);
volumePlot(row_indices, T, V0, policies_TS_OU, alpha, 'TS-OU', index_V0);

%% 5. OUT PRICE, CHECK FOR CONVERGENCE OF LS

%% Simulate the prices with a different seed
seedOut = 4;
% Simulation of the underlying for the OU case, low volatility.
XsL_OUT = spotSimulation('OU', paramsOUL, numberSimulations, 365, T, 0, seedOut, 0);
X_OUL_OUT = XsL_OUT(1:numberSimulations, :); XAV_OUL_OUT = XsL_OUT(numberSimulations+1:end, :);
S_OUL_OUT = S0*exp(X_OUL_OUT);
SAV_OUL_OUT = S0*exp(XAV_OUL_OUT);

% Simulation of the underlying for the OU case, high volatility.
XsH_OUT = spotSimulation('OU', paramsOUH, numberSimulations, 365, T, 0, seedOut, 0);
X_OUH_OUT = XsH_OUT(1:numberSimulations, :); XAV_OUH_OUT = XsH_OUT(numberSimulations+1:end, :);
S_OUH_OUT = S0*exp(X_OUH_OUT);
SAV_OUH_OUT = S0*exp(XAV_OUH_OUT);

% Simulation of the underlying for the OU-NTS symmetric case.
X_OU_NTSsymm_OUT = spotSimulation('OU-NTS', paramsNTSsymm, numberSimulations, 365, T, 16, seedOut, 1e-8);
S_OU_NTSsymm_OUT = S0*exp(X_OU_NTSsymm_OUT);

% Simulation of the underlying for the OU-NTS asymmetric case.
X_OU_NTSasymm_OUT = spotSimulation('OU-NTS', paramsNTSasymm, numberSimulations, 365, T, 16, seedOut, 1e-8);
S_OU_NTSasymm_OUT = S0*exp(X_OU_NTSasymm_OUT);

% Simulation of the underlying for the NTS-OU symmetric case.
X_NTS_OUsymm_OUT = spotSimulation('NTS-OU', paramsNTSsymm, numberSimulations, 365, T, 24, seedOut, 1e-10);
S_NTS_OUsymm_OUT = S0*exp(X_NTS_OUsymm_OUT);

% Simulation of the underlying for the NTS-OU asymmetric case.
X_NTS_OUasymm_OUT = spotSimulation('NTS-OU', paramsNTSasymm, numberSimulations, 365, T, 24, seedOut, 1e-10);
S_NTS_OUasymm_OUT = S0*exp(X_NTS_OUasymm_OUT);

% Simulation of the underlying for the OU-TS case.
X_OU_TS_OUT = spotSimulation('OU-TS', paramsTS, numberSimulations, 365, T, 18, seedOut, 1e-10);
S_OU_TS_OUT = S0*exp(X_OU_TS_OUT);

% Simulation of the underlying for the TS-OU case.
X_TS_OU_OUT = spotSimulation('TS-OU', paramsTS, numberSimulations, 365, T, 25, seedOut, 1e-12);
S_TS_OU_OUT = S0*exp(X_TS_OU_OUT);

%% Compute the final cashflows for the backward process by using out simulations.
% OU, low and high volatilities, normal and AV.
cashflows_OUL_OUT_T = penFunc(S_OUL_OUT(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVL_OUT_T = penFunc(SAV_OUL_OUT(:,end), ones(numberSimulations,1)*dV');
cashflows_OUH_OUT_T = penFunc(S_OUH_OUT(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AVH_OUT_T = penFunc(SAV_OUH_OUT(:,end), ones(numberSimulations,1)*dV');

% OU-NTS, symmetric and asymmetric.
cashflows_OU_NTSsymm_OUT_T = penFunc(S_OU_NTSsymm_OUT(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_NTSasymm_OUT_T = penFunc(S_OU_NTSasymm_OUT(:,end), ones(numberSimulations,1)*dV');

% NTS-OU symmetric and asymmetric.
cashflows_NTS_OUsymm_OUT_T = penFunc(S_NTS_OUsymm_OUT(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUasymm_OUT_T = penFunc(S_NTS_OUasymm_OUT(:,end), ones(numberSimulations,1)*dV');

% OU-TS.
cashflows_OU_TS_OUT_T = penFunc(S_OU_TS_OUT(:,end), ones(numberSimulations,1)*dV');

% TS-OU.
cashflows_TS_OU_OUT_T = penFunc(S_TS_OU_OUT(:,end), ones(numberSimulations,1)*dV');

%% Apply priceOut
% OU low volatility.
cashflows_OUL_OUT = priceOut(S_OUL_OUT, cashflows_OUL_OUT_T, regressors_OUL, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVL_OUT = priceOut(SAV_OUL_OUT, cashflows_OU_AVL_OUT_T, regressors_OU_AVL, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU high volatility.
cashflows_OUH_OUT = priceOut(S_OUH_OUT, cashflows_OUH_OUT_T, regressors_OUH, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVH_OUT = priceOut(SAV_OUH_OUT, cashflows_OU_AVH_OUT_T, regressors_OU_AVH, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-NTS, symmetric and asymmetric.
cashflows_OU_NTSsymm_OUT = priceOut(S_OU_NTSsymm_OUT, cashflows_OU_NTSsymm_OUT_T, regressors_OU_NTSsymm, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_NTSasymm_OUT = priceOut(S_OU_NTSasymm_OUT, cashflows_OU_NTSasymm_OUT_T, regressors_OU_NTSasymm, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% NTS-OU, symmetric and asymmetric.
cashflows_NTS_OUsymm_OUT = priceOut(S_NTS_OUsymm_OUT, cashflows_NTS_OUsymm_OUT_T, regressors_NTS_OUsymm, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OUasymm_OUT = priceOut(S_NTS_OUasymm_OUT, cashflows_NTS_OUasymm_OUT_T, regressors_NTS_OUasymm, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-TS.
cashflows_OU_TS_OUT = priceOut(S_OU_TS_OUT, cashflows_OU_TS_OUT_T, regressors_OU_TS, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% TS-OU.
cashflows_TS_OU_OUT = priceOut(S_TS_OU_OUT, cashflows_TS_OU_OUT_T, regressors_TS_OU, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

%% Compute the mean, the ci and the std error with the out prices

% OU low volatility.
[price_OUL_OUT, ~, price_OUL_OUT_CI] = normfit(0.5*(cashflows_OUL_OUT(:,index_V0) + cashflows_OU_AVL_OUT(:,index_V0)));
OUT_OUL_STD = standardError(0.5*(cashflows_OUL_OUT(:,index_V0) + cashflows_OU_AVL_OUT(:,index_V0)));

% OU high volatility.
[price_OUH_OUT, ~, price_OUH_OUT_CI] = normfit(0.5*(cashflows_OUH_OUT(:,index_V0) + cashflows_OU_AVH_OUT(:,index_V0)));
OUT_OUH_STD = standardError(0.5*(cashflows_OUH_OUT(:,index_V0) + cashflows_OU_AVH_OUT(:,index_V0)));

% OU-NTS symmetric.
[price_OU_NTSsymm_OUT, ~, price_OU_NTSsymm_OUT_CI] = normfit(cashflows_OU_NTSsymm_OUT(:, index_V0));
OUT_OU_NTSsymm_STD = standardError(cashflows_OU_NTSsymm_OUT(:, index_V0));

% OU-NTS asymmetric.
[price_OU_NTSasymm_OUT, ~, price_OU_NTSasymm_OUT_CI] = normfit(cashflows_OU_NTSasymm_OUT(:, index_V0));
OUT_OU_NTSasymm_STD = standardError(cashflows_OU_NTSasymm_OUT(:, index_V0));

% NTS-OU symmetric.
[price_NTS_OUsymm_OUT, ~, price_NTS_OUsymm_OUT_CI] = normfit(cashflows_NTS_OUsymm_OUT(:, index_V0));
OUT_NTS_OUsymm_STD = standardError(cashflows_NTS_OUsymm_OUT(:, index_V0));

% NTS-OU asymmetric.
[price_NTS_OUasymm_OUT, ~, price_NTS_OUasymm_OUT_CI] = normfit(cashflows_NTS_OUasymm_OUT(:, index_V0));
OUT_NTS_OUasymm_STD = standardError(cashflows_NTS_OUasymm_OUT(:, index_V0));

% OU-TS.
[price_OU_TS_OUT, ~, price_OU_TS_OUT_CI] = normfit(cashflows_OU_TS_OUT(:, index_V0));
OUT_OU_TS_STD = standardError(cashflows_OU_TS_OUT(:, index_V0));

% TS-OU.
[price_TS_OU_OUT, ~, price_TS_OU_OUT_CI] = normfit(cashflows_TS_OU_OUT(:, index_V0));
OUT_TS_OU_STD = standardError(cashflows_TS_OU_OUT(:, index_V0));