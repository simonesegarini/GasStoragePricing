% Script used to generates plots to show the impact of different sets of
% parameters on the long mean term of the OU process.

clear all, close all, clc

numberSimulations = 500;
T = 365;
S0 = 14.88;
seed = 2;
paramsOUL = [0.0315, 0.05, -0.3, 0.1;
             0.0315, 0.05, 0.3, 0.1;
             0.0315, 0.05, 0.3, -0.1;
             0.0315, 0.05, -0.3, -0.1]; % Parameters for the OU case, low volatility.

% Number of different parameter sets.
numParams = size(paramsOUL, 1);

% Plot trajectories for each parameter set.
figure;
hold on;

% Iterate over each parameter set.
for i = 1:numParams
    % Extract parameters for the ith set.
    params = paramsOUL(i, :);
    
    % Simulate the spot process.
    XsL = spotSimulation('OU', params, numberSimulations, 365, T, 0, seed, 0);
    
    % Take onyl the normal OU, not the AV part.
    X_OUL = XsL(1:numberSimulations, :);
    S_OUL = S0 * exp(X_OUL);
    
    % Plot a single trajectory for the ith parameter set.
    plot(S_OUL(1, :), 'DisplayName', sprintf('Parameters: A = %.1f, B = %.1f', params(3), params(4)));
end

% Details for the plot.
xlabel('Days')
ylabel('Spot Price')
title('Spot Price Trajectories for Different Parameters')
legend('show', 'Location','northwest')
grid on
hold off