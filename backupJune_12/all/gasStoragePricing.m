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

% OU
paramsOU = [0.0315, 0.05, 0]; % Parameters for the OU case
% params = @(t) [0.0945, 0.05, 0*t];

% OU-NTS
paramsNTS = [0.2, 0.2162, 0.201, 0.256, 0]; % Params for the OU-NTS case

% OU-TS
paramsTS = [0.4, 0.1, 2.5, 3.5, 0.5, 1, 0]; % Params for the OU-TS case

% Number of simulations
M = 500; 

% Discretization for the backward induction
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Simulation of the underlying for the OU case
Xs = spotSimulation('OU', paramsOU, M, 365, T, 0, 1);
X_OU = Xs(1:M, :); XAV_OU = Xs(M+1:end, :);
S_OU = S0*exp(X_OU);
SAV_OU = S0*exp(XAV_OU);

% Simulation of the underlying for the OU-NTS case
tic
X_OU_NTS = spotSimulation('OU-NTS', paramsNTS, M, 365, T, 16, 1);
timeOU_NTS = toc;
S_OU_NTS = S0*exp(X_OU_NTS);

% Simulation of the underlying for the OU-NTS case
tic
X_NTS_OU = spotSimulation('NTS-OU', paramsNTS, M, 365, T, 25, 1);
timeNTS_OU = toc;
S_NTS_OU = S0*exp(X_NTS_OU);

% Simulation of the underlying for the OU-TS case
tic
X_OU_TS = spotSimulation('OU-TS', paramsTS, M, 365, T, 21, 1);
timeOU_TS = toc;
S_OU_TS = S0*exp(X_OU_TS);

%% Plot to see some paths
simulations = [1, 5, 25, 50];
plotSpot(simulations, S_OU, 'OU')
plotSpot(simulations, SAV_OU, 'OU AV')
plotSpot(simulations, S_OU_NTS, 'OU-NTS')
plotSpot(simulations, S_NTS_OU, 'NTS-OU')
plotSpot(simulations, S_OU_TS, 'OU-TS')

%% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
% Matrices with the cashflows
cashflows_OU = penFunc(S_OU(:,end), ones(M,1)*dV');
cashflows_OU_AV = penFunc(SAV_OU(:,end), ones(M,1)*dV');
cashflows_OU_NTS = penFunc(S_OU_NTS(:,end), ones(M,1)*dV');
cashflows_NTS_OU = penFunc(S_NTS_OU(:,end), ones(M,1)*dV');
cashflows_OU_TS = penFunc(S_OU_TS(:,end), ones(M,1)*dV');

% Matrices for withdraw/injection limits
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

%% 3. APPLY BACKWARD INDUCTION
cashflows_OU = priceIn(S_OU, cashflows_OU, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
cashflows_OU_AV = priceIn(SAV_OU, cashflows_OU_AV, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);

cashflows_OU_NTS = priceIn(S_OU_NTS, cashflows_OU_NTS, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
cashflows_NTS_OU = priceIn(S_NTS_OU, cashflows_NTS_OU, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
cashflows_OU_TS = priceIn(S_OU_TS, cashflows_OU_TS, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows, policies, regression] = backwardInduction(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_AV, policies_AV, regressionAV] = backwardInduction(SAV, cashflows_AV, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);

% NEED POLICIES FOR THIS, RUN BACKWARDINDUCTION NOT PRICEIN!
% %% PLOT OF THE VOLUME EVOLUTION
% volumePlot(row_indices, T, V0, policies, alpha, vol, index_V0);

%% 4. PRICE
% Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
[price_OU_IN, IN_OU_STD, price_OU_IN_CI] = normfit(0.5*(cashflows_OU(:,index_V0) + cashflows_OU_AV(:,index_V0)))
[price_OU_NTS_IN, IN_OU_NTS_STD, price_OU_NTS_IN_CI] = normfit(cashflows_OU_NTS(:, index_V0))
[price_NTS_OU_IN, IN_NTS_OU_STD, price_NTS_OU_IN_CI] = normfit(cashflows_NTS_OU(:, index_V0))
[price_OU_TS_IN, IN_OU_TS_STD, price_OU_TS_IN_CI] = normfit(cashflows_OU_TS(:, index_V0))

% Plot the distribution of the IN price
figure;
histogram(0.5*(cashflows_OU(:,index_V0) + cashflows_OU_AV(:,index_V0)), 'NumBins', M/10)
title(['Value distribution (' num2str(M) ' paths), OU', ', sigma = ', num2str(paramsOU(1)*100), '%, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

figure;
histogram(cashflows_OU_NTS, 'NumBins', M/10)
title(['Value distribution (' num2str(M) ' paths), OU-NTS, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

figure;
histogram(cashflows_NTS_OU, 'NumBins', M/10)
title(['Value distribution (' num2str(M) ' paths), NTS-OU, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')

figure;
histogram(cashflows_OU_TS, 'NumBins', M/10)
title(['Value distribution (' num2str(M) ' paths), OU-TS, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')