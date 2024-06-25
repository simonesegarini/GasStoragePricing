function X = fgmc(model, params, N, M, T, Mfft, seed, toll)
% Finite General Monte Carlo algorithm for OU-Lévy and Lévy-OU processes
%
% INPUT:
% model:                char for model selection
% params:               vector with the parameters of the model
% N:                    number of simulations
% M:                    number of timesteps
% T:                    time horizon
% Mfft:                 FFT hyperparameter
% seed:                 set the seed of the simulation
% toll:                 tollerance for CDF selection
%
% OUTPUT:
% X:                    logprices

% Set given seed.
rng(seed)

% Parameters for the discretization.
dt = T/M;
alpha = params(1); b = params(2);

U = rand(N, M);
X = zeros(N,M+1);

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
for j = 1:M
    % Compute Zt based on the type of process.
    if strcmp(activity, 'Infinite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcIA(U(:,j), params, Mfft, dt, model, activity, toll);
        
    elseif strcmp(activity, 'Finite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcFA(U(:,j), params, Mfft, dt, model, activity, toll);
    end
end

end