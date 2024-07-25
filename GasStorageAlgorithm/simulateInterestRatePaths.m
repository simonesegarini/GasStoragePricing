function r = simulateInterestRatePaths(model, r0, paramsR, nSim, nSteps, T, seed)
% Simulate interest rate paths using Vasicek or CIR model.
%
% INPUT:
% model:                char for model selection
% paramsR:              vector with the parameters of the model
% r0:                   initial value for the IR
% nSim:                 number of simulations
% nSteps:               number of timesteps
% T:                    time horizon
% alg:                  algorithm selection (optional)
%
% OUTPUT:
% r:                    matrix of simulated interest rate paths

% Set seed if given.
if nargin == 7
    rng(seed)
end

% Time increment.
dt = T / nSteps;

% Assign parameters.
kappa = paramsR(1);
theta = paramsR(2);
sigma = paramsR(3);

% Preallocate the array for interest rate paths
r = zeros(nSim, nSteps + 1);
r(:, 1) = r0;
Z = randn(nSim, nSteps+1);


% Forward simulation of the interest rates.
switch model
    case 'Vasicek'

        % Vasicek model: dr = kappa*(theta - r)*dt + sigma*dW.
        for j=1:nSteps
            r(:, j+1) = r(:, j) + kappa.*(theta - r(:,j)).*dt + sigma.*sqrt(dt).*Z(:, j);
        end
    case 'CIR'

        % CIR model: dr = kappa*(theta - r)*dt + sigma*sqrt(r)*dW.
        for j=1:nSteps
            r(:, j+1) = r(:, j) + kappa.*(theta - r(:,j)).*dt + sigma.*sqrt(dt).*sqrt(r(:, j)).*Z(:, j);
        end
    otherwise
            error('Invalid model selection. Choose either ''Vasicek'' or ''CIR''.');
end
end
