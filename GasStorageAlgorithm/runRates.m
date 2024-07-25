% Script used to validate the function to simulate rates and to get the
% plots used in the report.

clear all, close all, clc

nPaths = 500;
nSteps = 365;
T = 365; % Total time in years
r0 = 0.02; % Initial rate

% Parameters for Vasicek model
kappaV = 0.3/T; % Speed of mean reversion
thetaV = 0.05; % Long-term mean rate
sigmaV = 0.02/sqrt(T); % Volatility

% Parameters for CIR model
kappaC = 0.2/T; % Speed of mean reversion
thetaC = 0.04; % Long-term mean rate
sigmaC = 0.1/sqrt(T); % Volatility

% Simulate interest rate paths
rVasicek = simulateInterestRatePaths('Vasicek', r0, [kappaV, thetaV, sigmaV], nPaths, nSteps, T, 2);
rCIR = simulateInterestRatePaths('CIR', r0, [kappaC, thetaC, sigmaC], nPaths, nSteps, T, 2);

% Plot the first 10 paths
figure;
plot((0:nSteps)/nSteps, rVasicek(1:10, :));
title(['Interest Rate Paths using Vasicek Model']);
xlabel('Time (years)');
ylabel('Interest Rate');
grid on;

figure;
plot((0:nSteps)/nSteps, rCIR(1:10, :));
title(['Interest Rate Paths using CIR Model']);
xlabel('Time (years)');
ylabel('Interest Rate');
grid on