function Xs = fgmc(model, params, N, M, T, Mfft, seed, scaling_fact)
% Finite General Monte Carlo algorithm for OU-Lévy and Lévy-OU processes.
%
% INPUT:
% model:                char for model selection
% params:               vector with the parameters of the model
% N:                    number of simulations
% M:                    number of timesteps
% T:                    time horizon
% Mfft:                 FFT hyperparameter
% seed:                 set the seed of the simulation
% scaling_fact:         parameter to stretch the CDF and handle low alphas
%
% OUTPUT:
% X:                    logprices

% Set given seed.
rng(seed)

% Check if there is a given scaling factor to handle low alphas.
if nargin < 8
    scaling_fact = 1;
end

% Parameters for the discretization.
dt = T/M;
alpha = params(1); b = params(2);

U = rand(N, M);
X = zeros(N,M+1);
XAV = zeros(N, M+1);

% Check if the model has finite/infinite activity and finite/infinite
% variation.
switch model
    case 'OU-NTS'
        if alpha < 0
            activity = 'Finite';
        else
            activity = 'Infinite';
        end
    case 'OU-TS'
        if alpha < 0
            activity = 'Finite';
        else
            activity = 'Infinite';
        end
    case {'NTS-OU', 'TS-OU'}
        if alpha == 0
            activity = 'Finite';
        else
            activity = 'Infinite';
        end
end

% Iteration in discretized time to update the simulations.

% Create approximated CDF.
[xgrid_hat, CDF_hat] = createCDF(params, Mfft, dt, model, activity, scaling_fact);

% Compute Zt based on the type of process.
if strcmp(activity, 'Infinite')
    for j = 1:M
        % V1 di AV
%         Zt = fgmcIA(xgrid_hat, CDF_hat, U(:,j));
%         X(:, j+1) = exp(-b.*dt).*X(:, j) + Zt;
%         XAV(:, j+1) = exp(-b.*dt).*XAV(:, j) - Zt;

        % V2 di AV
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcIA(xgrid_hat, CDF_hat, U(:,j));
        XAV(:, j+1) = exp(-b.*dt).*XAV(:, j) + fgmcIA(xgrid_hat, CDF_hat, 1-U(:,j));
    end
elseif strcmp(activity, 'Finite')
    for j = 1:M
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcFA(xgrid_hat, CDF_hat, U(:,j), dt, params, model);
    end
end
% checkAntitheticCovariance(X, XAV);
Xs = [X; XAV];
end