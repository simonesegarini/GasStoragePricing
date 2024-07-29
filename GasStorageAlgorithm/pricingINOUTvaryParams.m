clear all, close all, clc
warning('off', 'all')

%% STANDARD STORAGE CONTRACT DATA
% Cost of injection parameters.
a1 = 0; b1 = 0;

% Cost of selling parameters.
a2 = 0; b2 = 0;


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

% Parameters.
% OU-TS.
paramsTS = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0];

% Number of simulations.
numberSimulations = [10000, 15000, 25000];

% Discretization for the backward induction.
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

% Choose the model.
model = 'OU-TS';

% Allocate the values to be saved.
price_IN = zeros(3,1); IN_STD = zeros(3,1); price_IN_CI = zeros(3,2);
price_OUT = zeros(3,1); OUT_STD = zeros(3,1); price_OUT_CI = zeros(3,2);


for z=1:numel(numberSimulations)

    numSims = numberSimulations(z);

    % Discount rate.
    delta = zeros(numSims, T+1);

    % Simulation of the underlying for the OU-NTS case
    disp(['Now doing underling simulations for #sim = ', num2str(numSims)])
    Xs = spotSimulation(model, paramsTS, numSims, 365, T, 16, 1, 1);
    X = Xs(1:numSims, :); %XAV = Xs_OU_TS(numberSimulations(simIt)+1:end, :);
    S = S0*exp(X);
    %SAV = S0*exp(XAV);
    clear Xs;

    % 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
    % Matrices with the cashflows
    cashflows_IN = penFunc(S(:,end), ones(numSims,1)*dV');
    
    disp(['Now applying backward induction for #sim = ', num2str(numSims)])
    % 3. APPLY BACKWARD INDUCTION
    [cashflows_IN, ~, regressors_IN] = backwardInduction(S, cashflows_IN, payoff, N, numSims, delta, alpha, T, maxInjection, maxWithdraw, 5, 'polynomial');
    
    disp(['Now princing (in) for #sim = ', num2str(numSims)])
    % 4. PRICE IN
    [price_IN(z,:), ~ , price_IN_CI(z,:)] = normfit(cashflows_IN(:, index_V0));
    IN_STD(z,:) = standardError(cashflows_IN(:, index_V0));
    clear cashflows_IN;

    % 5. PRICE OUT
    disp(['Now doing out procedure for simulations for #sim = ', num2str(numSims)])
    Xs_OUT = spotSimulation(model, paramsTS, numSims, 365, T, 16, 2, 1);
    X_OUT = Xs_OUT(1:numSims, :);
    S_OUT = S0*exp(X_OUT);
    clear Xs_OUT;

    cashflows_OUT = penFunc(S_OUT(:,end), ones(numSims,1)*dV');
    cashflows_OUT = priceOut(S_OUT, cashflows_OUT, regressors_IN, payoff, N, numSims, delta, alpha, T, maxInjection, maxWithdraw, 5, 'polynomial');
    
    [price_OUT(z,:), ~ , price_OUT_CI(z,:)] = normfit(cashflows_OUT(:, index_V0));
    OUT_STD(z,:) = standardError(cashflows_OUT(:, index_V0));
    clear cashflows_OUT;

end