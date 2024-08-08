clear all, close all, clc

%% Set up FLAGS for processes.
OU = 1; 
OU_TS_FGMC = 21; TS_OU_FGMC = 31; 
OU_TS_SSR = 22; TS_OU_SR = 32;
OU_TS_DR = 23; TS_OU_DR = 33;

OU_NTS = 4; NTS_OU = 5;

%% Set up FLAGS for regression.
POLYNOMIAL = 1;
BSPLINE = 2;

%% Set up FLAGS for rates simluation.
NULL = 0;
CIR = 1;
VASICEK = 2;

%% Define the struct for storage pricing.

contract.SpotInitial = 14.88;
contract.Simulations = 500; % Number of paths simulated.
contract.Maturity = 365; % Maturity expressed in days.
contract.Steps = 365; % Steps for simulation, is equal to Maturity for daily simulations.
contract.Seed = 2; % Seed for the rng, NOT mandatory.

contract.Volume.Vmin = 0; % Minimum volume of the storage.
contract.Volume.Vmax = 250000; % Maximum volume of the storage.
contract.Volume.V0 = 100000; % Volume of the storage when starting simulations.
contract.Volume.VT = 100000; % Volume of the storage needed at maturity.
contract.Volume.MaxWithdraw = -7500; % Max allowed withdraw.
contract.Volume.MaxInjection = 2500; % Max allowed injection.
contract.Volume.Discretization = 2500; % Discretization of the volume for regression.

contract.Costs.TransactionCost = 0; %a1
contract.Costs.SpreadCost = 0; %b1
contract.Costs.TransactionProfit = 0; %a2
contract.Costs.SpreadProfit = 0; %b2

% SPOT PRICES SIMULATIONS.
% Process available: OU, OU-NTS, NTS-OU, OU-TS, TS-OU with also all special
% cases.
% If OU is selected, parameters should be given as [sigma, kappa, A, B] as
% the mean is modelled as A.*cos(2.*pi.*t./T) + B.
% If OU-NTS or NTS-OU are selected, parameters should be given as 
% [alpha, b, sigma, k, theta].
% If OU-TS or TS-OU are selected, parameters should be given as 
% [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c].
contract.GasProcess = OU_TS_SSR;
contract.GasParameters = [0.7, 0.1, 2.5, 3.5, 0.5, 1, 0];

% RATES SIMULATIONS.
% CIR, Vasicek, NULL (zeros) simulations available, give proper parameters.
% If NULL selected then initial rate and parameters are ignored and can 
% give 0s as input, otherwise give parameters in the format 
% [kappa, theta, sigma] for CIR and Vasicek rates model.
contract.Rates.InitialRate = 0;
contract.Rates.Process = NULL; 
contract.Rates.Parameters = [0, 0, 0]; 

% REGRESSION SETTINGS.
% Regressions available: polynomial and b-spline.
contract.Regression = POLYNOMIAL;
contract.Order = 4;

%% PRICE
[contract.Price, contract.CI, contract.SE, contract.ExecutionTime] = storagePricing(contract);

disp(['Contract price: ', num2str(contract.Price / 1e6, '%.3f'), ' million euros.'])
disp(['Contract price CI: (', num2str(contract.CI(1) / 1e6, '%.3f'), ', ', num2str(contract.CI(2) / 1e6, '%.3f'), ') million euros.'])
disp(['Contract price SE: ', num2str(contract.SE / 1e3, '%.3f'), ' thousand euros.'])
disp(['Running time: ', num2str(contract.ExecutionTime, '%.2f'), ' seconds.'])