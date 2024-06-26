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
paramsOU = [0.0315, 0.05, 0]; % Parameters for the OU case

% OU-NTS.
paramsNTS = [0.6, 0.2162, 0.201, 0.256, 0]; % Params for the OU-NTS case

% OU-TS.
paramsTS = [0.8, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case

% Number of simulations.
numberSimulations = 500; 

% Discretization for the backward induction.
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Simulation of the underlying for the OU case.
Xs = spotSimulation('OU', paramsOU, numberSimulations, 365, T, 0, seed, 0);
X_OU = Xs(1:numberSimulations, :); XAV_OU = Xs(numberSimulations+1:end, :);
S_OU = S0*exp(X_OU);
SAV_OU = S0*exp(XAV_OU);

% Simulation of the underlying for the OU-NTS case.
tic
X_OU_NTS = spotSimulation('OU-NTS', paramsNTS, numberSimulations, 365, T, 16, seed, 1e-8);
timeOU_NTS = toc;
S_OU_NTS = S0*exp(X_OU_NTS);

% Simulation of the underlying for the NTS-OU case.
tic
X_NTS_OU = spotSimulation('NTS-OU', paramsNTS, numberSimulations, 365, T, 24, seed, 1e-10);
timeNTS_OU = toc;
S_NTS_OU = S0*exp(X_NTS_OU);

% Simulation of the underlying for the OU-TS case.
tic
X_OU_TS = spotSimulation('OU-TS', paramsTS, numberSimulations, 365, T, 18, seed, 1e-10);
timeOU_TS = toc;
S_OU_TS = S0*exp(X_OU_TS);

% Simulation of the underlying for the TS-OU case.
tic
X_TS_OU = spotSimulation('TS-OU', paramsTS, numberSimulations, 365, T, 25, seed, 1e-12);
timeTS_OU = toc;
S_TS_OU = S0*exp(X_TS_OU);

%% Plot to see some paths
simulations = [1, 5, 25, 50];
plotSpot(simulations, S_OU, 'OU')
plotSpot(simulations, SAV_OU, 'OU AV')
plotSpot(simulations, S_OU_NTS, 'OU-NTS')
plotSpot(simulations, S_NTS_OU, 'NTS-OU')
plotSpot(simulations, S_OU_TS, 'OU-TS')
plotSpot(simulations, S_TS_OU, 'TS-OU')

%% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
% Assign the final penalization cost to the cashflows matrices, they have a
% T in the end to refer as the final time step, created this variable to
% avoid new initialization each time we want to try a different approach.
cashflows_OU_T = penFunc(S_OU(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_AV_T = penFunc(SAV_OU(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_NTS_T = penFunc(S_OU_NTS(:,end), ones(numberSimulations,1)*dV');
cashflows_NTS_OU_T = penFunc(S_NTS_OU(:,end), ones(numberSimulations,1)*dV');
cashflows_OU_TS_T = penFunc(S_OU_TS(:,end), ones(numberSimulations,1)*dV');
cashflows_TS_OU_T = penFunc(S_TS_OU(:,end), ones(numberSimulations,1)*dV');

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

%% 3. APPLY BACKWARD INDUCTION
% Use IN-pricing to avoid using too much space for the matrices of policies
% and for the regressors, they are needed just for plot and OUT-price
% purpose.
cashflows_OUP4 = priceIn(S_OU, cashflows_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_AVP4 = priceIn(SAV_OU, cashflows_OU_AV_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

cashflows_OU_NTS = priceIn(S_OU_NTS, cashflows_OU_NTS_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_NTS_OU = priceIn(S_NTS_OU, cashflows_NTS_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_OU_TS = priceIn(S_OU_TS, cashflows_OU_TS_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');
cashflows_TS_OU = priceIn(S_TS_OU, cashflows_TS_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 4, 'polynomial');

%% 4. PRICE
% Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
[price_OU_IN_P4, IN_OU_STD_P4, price_OU_IN_CI_P4] = normfit(0.5*(cashflows_OUP4(:,index_V0) + cashflows_OU_AVP4(:,index_V0)))
[price_OU_NTS_IN, IN_OU_NTS_STD, price_OU_NTS_IN_CI] = normfit(cashflows_OU_NTS(:, index_V0))
[price_NTS_OU_IN, IN_NTS_OU_STD, price_NTS_OU_IN_CI] = normfit(cashflows_NTS_OU(:, index_V0))
[price_OU_TS_IN, IN_OU_TS_STD, price_OU_TS_IN_CI] = normfit(cashflows_OU_TS(:, index_V0))
[price_TS_OU_IN, IN_TS_OU_STD, price_TS_OU_IN_CI] = normfit(cashflows_TS_OU(:, index_V0))

% Plot the distribution of the IN price
% OU case
figure;
histogram(0.5*(cashflows_OUP4(:,index_V0) + cashflows_OU_AVP4(:,index_V0)), 'NumBins', numberSimulations/10)
title(['Value distribution (' num2str(numberSimulations) ' paths), OU', ', sigma = ', num2str(paramsOU(1)*100), '%, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

% OU-NTS case
figure;
histogram(cashflows_OU_NTS, 'NumBins', numberSimulations/10)
title(['Value distribution (' num2str(numberSimulations) ' paths), OU-NTS, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

% NTS-OU case
figure;
histogram(cashflows_NTS_OU, 'NumBins', numberSimulations/10)
title(['Value distribution (' num2str(numberSimulations) ' paths), NTS-OU, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

% OU-TS case
figure;
histogram(cashflows_OU_TS, 'NumBins', numberSimulations/10)
title(['Value distribution (' num2str(numberSimulations) ' paths), OU-TS, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

% TS-OU case
figure;
histogram(cashflows_TS_OU, 'NumBins', numberSimulations/10)
title(['Value distribution (' num2str(numberSimulations) ' paths), TS-OU, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

%% 5. DIFFERENT BASIS
cashflows_OUB5 = priceIn(S_OU, cashflows_OU_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 5, 'bspline');
cashflows_OU_AVB5 = priceIn(SAV_OU, cashflows_OU_AV_T, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw, 5, 'bspline');

[price_OU_IN_B5, IN_OU_STD_B5, price_OU_IN_CI_B5] = normfit(0.5*(cashflows_OUB5(:,index_V0) + cashflows_OU_AVB5(:,index_V0)))

%% 6. OUT PRICE, CHECK FOR CONVERGENCE OF LS
% OUT-price, uncomment to run the backwardInduction algorithm that gives
% the regressors and the policies, run the part of priceOut after
% simulating the prices with another seed.
%
% [cashflows, policies, regression] = backwardInduction(S, cashflows, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_AV, policies_AV, regressionAV] = backwardInduction(SAV, cashflows_AV, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_OU_NTS, policies_OU_NTS, regression_OU_NTS] = backwardInduction(S_OU_NTS, cashflows_OU_NTS, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_NTS_OU, policies_NTS_OU, regression_NTS_OU] = backwardInduction(S_NTS_OU, cashflows_NTS_OU, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_OU_TS, policies_OU_TS, regression_OU_TS] = backwardInduction(S_OU_TS, cashflows_OU_TS, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_TS_OU, policies_TS_OU, regression_TS_OU] = backwardInduction(S_TS_OU, cashflows_TS_OU, payoff, N, numberSimulations, delta, alpha, T, maxInjection, maxWithdraw);


% NEED POLICIES FOR THIS, RUN BACKWARDINDUCTION NOT PRICEIN!
% %% PLOT OF THE VOLUME EVOLUTION
% volumePlot(row_indices, T, V0, policies, alpha, vol, index_V0);