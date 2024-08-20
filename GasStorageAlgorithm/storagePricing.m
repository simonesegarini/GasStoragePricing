function [price, CI, SE, time] = storagePricing(contract)
% Simulate rates and underlying and then proceed to price a storage
% contract with the details contained in the input struct.
%
% INPUT:
% contract:             struct with the contract
%
% OUTPUT:
% price:                price of the contract
% CI:                   confidence interval of the contract
% SE:                   standard error of the contract
% ratesSimulated:       matrix of simulated rates
% spotSimulated:        matrix with the simulated underlyings

% Check if a seed is given.
if isfield(contract, 'Seed') && ~isempty(contract.Seed)
    seed = contract.Seed;
else
    seed = randi(10000); % Generate a random seed between 1 and 10000.
end

% Check general struct of the contract.
checkContract(contract)

% Assign cost of injection parameters and cost of selling parameters.
a1 = contract.Costs.TransactionCost; b1 = contract.Costs.SpreadCost;
a2 = contract.Costs.TransactionProfit; b2 = contract.Costs.SpreadProfit;

% Define payoff for each simulation.
payoff = @(s, deltaV) -((1+a1).*s + b1).*deltaV.*(deltaV > 0) ...
                - ((1-a2).*s - b2).*deltaV.*(deltaV < 0);

% Assign volume constraints and values.
Vmin = contract.Volume.Vmin;
Vmax = contract.Volume.Vmax;
V0 = contract.Volume.V0;
VT = contract.Volume.VT;

% Assign maximum and minimum allowed injection and withdrawal of gas in time.
Imin = @(v) max(Vmin - v, contract.Volume.MaxWithdraw);
Imax = @(v) min(Vmax - v, contract.Volume.MaxInjection);

% Definition of penalization function in order to force the final desired
% storage volume.
penFunc = @(s, v) -s.*abs(v-VT).^2;

% Assign discretization parameters for the simulations of S and r.
numberSimulations = contract.Simulations; 
T = contract.Maturity;
M = contract.Steps;

% Discretization for the backward induction.
alpha = contract.Volume.Discretization; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

% Matrices for withdraw/injection limits.
maxWithdraw = Imin(dV*ones(1,T));
maxInjection = Imax(dV*ones(1,T));

% SIMULATION OF RATES.
r0 = contract.Rates.InitialRate;
ratesModel = contract.Rates.Process; 
ratesParameters = contract.Rates.Parameters; 
ratesSimulated = simulateInterestRatePaths(ratesModel, r0, ratesParameters, numberSimulations, M, T, seed);

% Check if the process allows AV method.
gasProcess = contract.GasProcess;
parameters = contract.GasParameters;
AV = logical(contract.AV);
if (gasProcess == 22 || gasProcess == 23 || gasProcess == 32 || gasProcess == 33)
    AV = false;
end

tic
if parameters(1) >= 1 && (gasProcess == 22 || gasProcess == 23 || gasProcess == 32 || gasProcess == 33)
    error('DR & SSR are not available for values of alpha greater or equal than 1.')
else
    if AV
        % SIMULATION OF UNDERLYINGS.
        [spotSimulated, spotSimulatedAV] = spotSimulation(contract.SpotInitial, gasProcess, parameters, numberSimulations, M, T, 16, contract.GPU, contract.STRETCH, seed);

        % BACKWARD INDUCTION TO PRICE.
        cashflows_T = penFunc(spotSimulated(:,end), ones(numberSimulations,1)*dV');
        cashflowsAV_T = penFunc(spotSimulatedAV(:,end), ones(numberSimulations,1)*dV');

        switch contract.Regression
            case 1
                cashflows = priceIn(spotSimulated, cashflows_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'polynomial');
                cashflowsAV = priceIn(spotSimulatedAV, cashflowsAV_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'polynomial');
            case 2
                cashflows = priceIn(spotSimulated, cashflows_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'bspline');
                cashflowsAV = priceIn(spotSimulatedAV, cashflowsAV_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'bspline');
            otherwise
                error('Select a valid regression for backward induction.')
        end

        [price, ~, CI] = normfit(0.5*(cashflows(:,index_V0) + cashflowsAV(:,index_V0)));
        SE = standardError(0.5*(cashflows(:,index_V0) + cashflowsAV(:,index_V0)));


    elseif ~AV
        % SIMULATION OF UNDERLYINGS.
        spotSimulated = spotSimulation(contract.SpotInitial, gasProcess, parameters, numberSimulations, M, T, 16, contract.GPU, seed);

        % BACKWARD INDUCTION TO PRICE.
        cashflows_T = penFunc(spotSimulated(:,end), ones(numberSimulations,1)*dV');

        switch contract.Regression
            case 1
                cashflows = priceIn(spotSimulated, cashflows_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'polynomial');
            case 2
                cashflows = priceIn(spotSimulated, cashflows_T, payoff, N, numberSimulations, ratesSimulated, alpha, T, maxInjection, maxWithdraw, contract.Order, 'bspline');
            otherwise
                error('Select a valid regression for backward induction.')
        end

        [price, ~, CI] = normfit(cashflows(:,index_V0));
        SE = standardError(cashflows(:,index_V0));
    end
end
time = toc;

if numberSimulations < 50
    simulationsPlot = [1, 2, 3, 4];
else
    simulationsPlot = [1, 5, 25, 50];
end

plotSpot(simulationsPlot, spotSimulated, gasProcess, 0)
if AV
    plotPriceDistribution(0.5*(cashflows(:,index_V0) + cashflowsAV(:,index_V0)), numberSimulations, gasProcess, parameters(1), 0)
elseif ~AV
    plotPriceDistribution(cashflows(:,index_V0), numberSimulations, gasProcess, parameters(1), 0)
end