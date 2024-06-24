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
k = 0.05; % MR parameter
mu = @(t) 0; % deterministic function of time (can be adjusted)
%mu = @(t) 0.05*t/365; 
sigma = [0.0315, 0.0945]; % low and high volatilities, MR parameter
M = 500; % number of simulations
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';

index_V0 = find(dV == V0);
vol = sigma(1);
[X, XAV] = spotSimulation(mu, vol, k, M, T+1, 1);
S = S0*exp(X);
SAV = S0*exp(XAV);

% %% PLOT OF THE SPOT EVOLUTION
% row_indices = [1, 10, 25, 50];
% plotSpot(row_indices, S, 'MC', vol);
% plotSpot(row_indices, SAV, 'MCAV', vol);

%% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
% Matrices with the cashflows
cashflows = penFunc(S(:,end), ones(M,1)*dV');
cashflows_AV = penFunc(SAV(:,end), ones(M,1)*dV');

% Matrices for withdraw/injection limits
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

%% 3. APPLY BACKWARD INDUCTION
tic
cashflows = priceIn(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
cashflows_AV = priceIn(SAV, cashflows_AV, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
time = toc;
% [cashflows, policies, regression] = backwardInduction(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
% [cashflows_AV, policies_AV, regressionAV] = backwardInduction(SAV, cashflows_AV, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);

% NEED POLICIES FOR THIS, RUN BACKWARDINDUCTION NOT PRICEIN!
% %% PLOT OF THE VOLUME EVOLUTION
% volumePlot(row_indices, T, V0, policies, alpha, vol, index_V0);

%% 4. PRICE
% Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
[price_IN, IN_STD, price_IN_CI] = normfit(0.5*(cashflows(:,index_V0) + cashflows_AV(:,index_V0)))

% Plot the distribution of the IN price
figure;
histogram(0.5*(cashflows(:,index_V0) + cashflows_AV(:,index_V0)), 'NumBins', M/10)
title(['Value distribution (' num2str(M) ' paths)', ', sigma = ', num2str(vol*100), '%, IN-SAMPLE']);
ylabel('Frequency')
xlabel('Value [Euro]')
