clear all, close all, clc

%% STANDARD STORAGE CONTRACT DATA
% Cost of injection parameters.
a1 = 0; b1 = 0;

% Cost of selling parameters.
a2 = 0; b2 = 0;

% Discount rate.
delta = 0.02;

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

% Discretizations for the simulations of S and r.
numberSimulations = 500; 
T = 365;
seed = 2;

%% 0. INTEREST RATES SIMULATION

% Initial rate.
r0 = 0.02;

% Parameters for Vasicek
kappaV = 0.3/T; % Speed of mean reversion
thetaV = 0.05; % Long-term mean rate
sigmaV = 0.02/sqrt(T); % Volatility
paramsRatesVasicek = [kappaV, thetaV, sigmaV];

% Parameters for CIR
kappaC = 0.2/T; % Speed of mean reversion
thetaC = 0.04; % Long-term mean rate
sigmaC = 0.1/sqrt(T); % Volatility
paramsRatesCIR = [kappaC, thetaC, sigmaC];

%rates = simulateInterestRatePaths('Vasicek', r0, paramsRatesVasicek, numberSimulations, 365, T, seed);
rates = simulateInterestRatePaths('CIR', r0, paramsRatesCIR, numberSimulations, 365, T, seed);
%rates = zeros(numberSimulations, T+1);

%% 1. SPOT PRICE SIMULATION 
S0 = 14.88;

% OU.
paramsOUL = [0.0315, 0.05, -0.3, 0.1]; % Parameters for the OU case, low volatility.
paramsOUH = [0.0945, 0.05, -0.3, 0.1]; % Parameters for the OU case, high volatility.

% OU-NTS.
paramsNTSsymm = [0.7, 0.2162, 0.201, 0.256, 0]; % Params for the OU-NTS symmetric case.
paramsNTSasymm = [0.7, 0.2162, 0.201, 0.256, 0.1]; % Params for the OU-NTS asymmetric case.

% OU-TS.
paramsTS = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case

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
Xs_OU_NTSsymm = spotSimulation('OU-NTS', paramsNTSsymm, numberSimulations, 365, T, 16, seed, 1e-8);
X_OU_NTSsymm = Xs_OU_NTSsymm(1:numberSimulations, :); XAV_OU_NTSsymm = Xs_OU_NTSsymm(numberSimulations+1:end, :);
timeOU_NTSsymm = toc;
S_OU_NTSsymm = S0*exp(X_OU_NTSsymm);
SAV_OU_NTSsymm = S0*exp(XAV_OU_NTSsymm);

% Simulation of the underlying for the OU-NTS asymmetric case.
tic
Xs_OU_NTSasymm = spotSimulation('OU-NTS', paramsNTSasymm, numberSimulations, 365, T, 16, seed, 1e-8);
X_OU_NTSasymm = Xs_OU_NTSasymm(1:numberSimulations, :); XAV_OU_NTSasymm = Xs_OU_NTSsymm(numberSimulations+1:end, :);
timeOU_NTSasymm = toc;
S_OU_NTSasymm = S0*exp(X_OU_NTSasymm);
SAV_OU_NTSasymm = S0*exp(XAV_OU_NTSasymm);

% Simulation of the underlying for the NTS-OU symmetric case.
tic
Xs_NTS_OUsymm = spotSimulation('NTS-OU', paramsNTSsymm, numberSimulations, 365, T, 24, seed, 1e-10);
X_NTS_OUsymm = Xs_NTS_OUsymm(1:numberSimulations, :); XAV_NTS_OUsymm = Xs_NTS_OUsymm(numberSimulations+1:end, :);
timeNTS_OUsymm = toc;
S_NTS_OUsymm = S0*exp(X_NTS_OUsymm);
SAV_NTS_OUsymm = S0*exp(XAV_NTS_OUsymm);

% Simulation of the underlying for the NTS-OU asymmetric case.
tic
Xs_NTS_OUasymm = spotSimulation('NTS-OU', paramsNTSasymm, numberSimulations, 365, T, 24, seed, 1e-10);
X_NTS_OUasymm = Xs_NTS_OUasymm(1:numberSimulations, :); XAV_NTS_OUasymm = Xs_NTS_OUasymm(numberSimulations+1:end, :);
timeNTS_OUasymm = toc;
S_NTS_OUasymm = S0*exp(X_NTS_OUasymm);
SAV_NTS_OUasymm = S0*exp(XAV_NTS_OUasymm);

% Simulation of the underlying for the OU-TS case.
tic
Xs_OU_TS = spotSimulation('OU-TS', paramsTS, numberSimulations, 365, T, 18, seed, 1e-10, 'fgmc');
X_OU_TS = Xs_OU_TS(1:numberSimulations, :); XAV_OU_TS = Xs_OU_TS(numberSimulations+1:end, :);
timeOU_TS = toc;
S_OU_TS = S0*exp(X_OU_TS);
SAV_OU_TS = S0*exp(XAV_OU_TS);

% Simulation of the underlying for the TS-OU case.
tic
Xs_TS_OU = spotSimulation('TS-OU', paramsTS, numberSimulations, 365, T, 25, seed, 1e-12, 'fgmc');
X_TS_OU = Xs_TS_OU(1:numberSimulations, :); XAV_TS_OU = Xs_TS_OU(numberSimulations+1:end, :);
timeTS_OU = toc;
S_TS_OU = S0*exp(X_TS_OU);
SAV_TS_OU = S0*exp(XAV_TS_OU);

%% Plot to see some paths
simulations = [1, 5, 25, 50];

% OU, low and high volatilities, normal and AV
plotSpot(simulations, S_OUL, 'OU-L', 0)
plotSpot(simulations, S_OUH, 'OU-H', 0)
plotSpot(simulations, SAV_OUL, 'OU-AVL', 0)
plotSpot(simulations, SAV_OUH, 'OU-AVH', 0)

% OU-NTS, symmetric and asymmetric
plotSpot(simulations, S_OU_NTSsymm, 'OU-NTSsymm', 0)
plotSpot(simulations, S_OU_NTSasymm, 'OU-NTSasymm', 0)
plotSpot(simulations, SAV_OU_NTSsymm, 'OU-NTSsymmAV', 0)
plotSpot(simulations, SAV_OU_NTSasymm, 'OU-NTSasymmAV', 0)

% NTS-OU symmetric and asymmetric
plotSpot(simulations, S_NTS_OUsymm, 'NTS-OUsymm', 0)
plotSpot(simulations, S_NTS_OUasymm, 'NTS-OUasymm', 0)
plotSpot(simulations, SAV_NTS_OUsymm, 'NTS-OUsymmAV', 0)
plotSpot(simulations, SAV_NTS_OUasymm, 'NTS-OUasymmAV', 0)

% OU-TS
plotSpot(simulations, S_OU_TS, 'OU-TS', 0)
plotSpot(simulations, SAV_OU_TS, 'OU-TSAV', 0)

% TS-OU
plotSpot(simulations, S_TS_OU, 'TS-OU', 0)
plotSpot(simulations, SAV_TS_OU, 'TS-OUAV', 0)

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
cashflows_OU_NTSsymmAV_T = penFunc(SAV_OU_NTSsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_NTSasymmAV_T = penFunc(SAV_OU_NTSasymm(:,end), ones(numberSimulations,1)*dV');

% NTS-OU symmetric and asymmetric
cashflows_NTS_OUsymm_T = penFunc(S_NTS_OUsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUasymm_T = penFunc(S_NTS_OUasymm(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUsymmAV_T = penFunc(SAV_NTS_OUsymm(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OUasymmAV_T = penFunc(SAV_NTS_OUasymm(:,end), ones(numberSimulations,1)*dV');

% OU-TS
cashflows_OU_TS_T = penFunc(S_OU_TS(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_TSAV_T = penFunc(SAV_OU_TS(:,end), ones(numberSimulations,1)*dV');

% TS-OU
cashflows_TS_OU_T = penFunc(S_TS_OU(:,end), ones(numberSimulations,1)*dV');
cashflows_TS_OUAV_T = penFunc(SAV_TS_OU(:,end), ones(numberSimulations,1)*dV');

%% 3. APPLY BACKWARD INDUCTION WITH PRICEIN
% Use IN-pricing to avoid using too much space for the matrices of policies
% and for the regressors, they are needed just for plot and OUT-price
% purpose.

% OU, low and high volatilities, normal and AV
cashflows_OUL = priceIn(S_OUL, cashflows_OUL_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVL = priceIn(SAV_OUL, cashflows_OU_AVL_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OUH = priceIn(S_OUH, cashflows_OUH_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVH = priceIn(SAV_OUH, cashflows_OU_AVH_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-NTS, symmetric and asymmetric
cashflows_OU_NTSsymm = priceIn(S_OU_NTSsymm, cashflows_OU_NTSsymm_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_NTSasymm = priceIn(S_OU_NTSasymm, cashflows_OU_NTSasymm_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_NTSsymmAV = priceIn(SAV_OU_NTSsymm, cashflows_OU_NTSsymmAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_NTSasymmAV = priceIn(SAV_OU_NTSasymm, cashflows_OU_NTSasymmAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% NTS-OU symmetric and asymmetric
cashflows_NTS_OUsymm = priceIn(S_NTS_OUsymm, cashflows_NTS_OUsymm_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OUasymm = priceIn(S_NTS_OUasymm, cashflows_NTS_OUasymm_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OUsymmAV = priceIn(SAV_NTS_OUsymm, cashflows_NTS_OUsymmAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OUasymmAV = priceIn(SAV_NTS_OUasymm, cashflows_NTS_OUasymmAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% OU-TS
cashflows_OU_TS = priceIn(S_OU_TS, cashflows_OU_TS_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_TSAV = priceIn(SAV_OU_TS, cashflows_OU_TSAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

% TS-OU
cashflows_TS_OU = priceIn(S_TS_OU, cashflows_TS_OU_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_TS_OUAV = priceIn(SAV_TS_OU, cashflows_TS_OUAV_T, payoff, N, numberSimulations, rates, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

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
[price_OU_NTSsymm_IN, ~, price_OU_NTSsymm_IN_CI] = normfit(0.5*(cashflows_OU_NTSsymm(:, index_V0) + cashflows_OU_NTSsymmAV(:, index_V0)));
IN_OU_NTSsymm_STD = standardError(0.5*(cashflows_OU_NTSsymm(:, index_V0) + cashflows_OU_NTSsymmAV(:, index_V0)));

% OU-NTS asymmetric
[price_OU_NTSasymm_IN, ~, price_OU_NTSasymm_IN_CI] = normfit(0.5*(cashflows_OU_NTSasymm(:, index_V0) + cashflows_OU_NTSasymmAV(:, index_V0)));
IN_OU_NTSasymm_STD = standardError(0.5*(cashflows_OU_NTSasymm(:, index_V0) + cashflows_OU_NTSasymmAV(:, index_V0)));

% NTS-OU symmetric
[price_NTS_OUsymm_IN, ~, price_NTS_OUsymm_IN_CI] = normfit(0.5*(cashflows_NTS_OUsymm(:, index_V0) + cashflows_NTS_OUsymmAV(:, index_V0)));
IN_NTS_OUsymm_STD = standardError(0.5*(cashflows_NTS_OUsymm(:, index_V0) + cashflows_NTS_OUsymmAV(:, index_V0)));

% NTS-OU asymmetric
[price_NTS_OUasymm_IN, ~, price_NTS_OUasymm_IN_CI] = normfit(0.5*(cashflows_NTS_OUasymm(:, index_V0) + cashflows_NTS_OUasymmAV(:, index_V0)));
IN_NTS_OUasymm_STD = standardError(0.5*(cashflows_NTS_OUasymm(:, index_V0) + cashflows_NTS_OUasymmAV(:, index_V0)));

% OU-TS
[price_OU_TS_IN, ~, price_OU_TS_IN_CI] = normfit(0.5*(cashflows_OU_TS(:, index_V0) + cashflows_OU_TSAV(:, index_V0)));
IN_OU_TS_STD = standardError(0.5*(cashflows_OU_TS(:, index_V0) + cashflows_OU_TSAV(:, index_V0)));

% TS-OU
[price_TS_OU_IN, ~, price_TS_OU_IN_CI] = normfit(0.5*(cashflows_TS_OU(:, index_V0) + cashflows_TS_OUAV(:, index_V0)));
IN_TS_OU_STD = standardError(0.5*(cashflows_TS_OU(:, index_V0) + cashflows_TS_OUAV(:, index_V0)));

%% Plots of distribution of storage prices
% Plot the distribution of the IN price (save the plots if the last inputs is 1)

% OU case, low volatility
plotPriceDistribution(0.5*(cashflows_OUL(:,index_V0) + cashflows_OU_AVL(:,index_V0)), numberSimulations, 'OU-L', paramsOUL(1)*100, 0);

% OU case, high volatility
plotPriceDistribution(0.5*(cashflows_OUH(:,index_V0) + cashflows_OU_AVH(:,index_V0)), numberSimulations, 'OU-H', paramsOUH(1)*100, 0);

% OU-NTS symmetric case
plotPriceDistribution(cashflows_OU_NTSsymm(:, index_V0), numberSimulations, 'OU-NTSsymm', paramsNTSsymm(1), 0);

% OU-NTS asymmetric case
plotPriceDistribution(cashflows_OU_NTSasymm(:, index_V0), numberSimulations, 'OU-NTSasymm', paramsNTSasymm(1), 0);

% NTS-OU symmetric case
plotPriceDistribution(cashflows_NTS_OUsymm(:, index_V0), numberSimulations, 'NTS-OUsymm', paramsNTSsymm(1), 0);

% NTS-OU asymmetric case
plotPriceDistribution(cashflows_NTS_OUasymm(:, index_V0), numberSimulations, 'NTS-OUasymm', paramsNTSasymm(1), 0);

% Note: for this cases it might be useful to remove the maximum value of
% the price vector in order to obtain graphs that are more readable.
% OU-TS case
plotPriceDistribution(cashflows_OU_TS(:, index_V0), numberSimulations, 'OU-TS', paramsTS(1), 0);

% TS-OU case
plotPriceDistribution(cashflows_TS_OU(:, index_V0), numberSimulations, 'TS-OU', paramsTS(1), 0);