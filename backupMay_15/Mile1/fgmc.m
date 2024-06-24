function [X, f] = fgmc(X0, params, N, M, T, model)

% DA IMPLEMENTARE POI:
% - CHECK CUMULANTS TS

% Parameters for the discretization
rng(2) % set seed for reproducibility
dt = T/M;
U = rand(N, M);
X = zeros(N,M+1); X(:,1) = X0;
f = zeros(1,M+1); f(:,1) = -LogCharFunc(-1i, 0, params, model);
Mfft = 16;
alpha = params(1);
b = params(2);

% Check if the model has finite/infinite activity and finite/infinite variation
if strcmp(model, 'NTS')
    if alpha < 0
        activity = 'Finite';
    elseif alpha >= 0 && alpha < 1
        activity = 'Infinite';
    end
elseif strcmp(model, 'TS')
    if alpha < 0
        activity = 'Finite';
    elseif alpha >= 0 && alpha < 2
        activity = 'Infinite';
    end
end

% Iteration in discretized time to update the simulations
for j = 1:M
    % Compute Zt based on the type of process
    if strcmp(activity, 'Infinite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcIA(U(:,j), params, Mfft, dt, model, activity);
    elseif strcmp(activity, 'Finite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcFA(U(:,j), params, Mfft, dt, model, activity);
    end
    f(:, j+1) = - LogCharFunc(-1i, dt*j, params, model); 
end

end