% Script used to validate the function to simulate rates and to get the plots used in the report.

clear all, close all, clc

% Parameters
nPaths = 500; % Number of paths to simulate
nSteps = 365; % Number of steps per path
T = 1; % Total time in years (adjusted to reflect actual time in years)
r0 = 0.02; % Initial interest rate

% Parameters for Vasicek model
kappaV = 0.3; % Speed of mean reversion
thetaV = 0.05; % Long-term mean rate
sigmaV = 0.02; % Volatility

% Parameters for CIR model
kappaC = 0.2; % Speed of mean reversion
thetaC = 0.04; % Long-term mean rate
sigmaC = 0.01; % Volatility

% Simulate interest rate paths using the Vasicek model
rVasicek = simulateInterestRatePaths('Vasicek', r0, [kappaV, thetaV, sigmaV], nPaths, nSteps, T, 2);

% Simulate interest rate paths using the CIR model
rCIR = simulateInterestRatePaths('CIR', r0, [kappaC, thetaC, sigmaC], nPaths, nSteps, T, 2);

% Define time vector for plotting
timeVector = linspace(0, T, nSteps + 1);

% Plot the first 10 paths for the Vasicek model
figure;
hold on;
plot(timeVector, rVasicek(1:10, :), 'LineWidth', 1.5);
title('Interest Rate Paths using Vasicek Model', 'FontSize', 14);
xlabel('Time (years)', 'FontSize', 12);
ylabel('Interest Rate', 'FontSize', 12);
grid on;
hold off;

% Plot the first 10 paths for the CIR model
figure;
hold on;
plot(timeVector, rCIR(1:10, :), 'LineWidth', 1.5);
title('Interest Rate Paths using CIR Model', 'FontSize', 14);
xlabel('Time (years)', 'FontSize', 12);
ylabel('Interest Rate', 'FontSize', 12);
grid on;
hold off;
