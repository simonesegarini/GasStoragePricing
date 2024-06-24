function X = fgmc(X0, params, N, M, T, model)

% Parameters for the discretization
rng(2) % set seed for reproducibility
dt = T/M;
U = rand(N, M);
X = zeros(N,M+1); X(:,1) = X0;
Mfft = 16;
alpha = params(1);
b = params(2);

% Check if the model has finite/infinite activity and finite/infinite variation
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

% Iteration in discretized time to update the simulations
for j = 1:M
    % Compute Zt based on the type of process
    if strcmp(activity, 'Infinite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcIA(U(:,j), params, Mfft, dt, model, activity);
    elseif strcmp(activity, 'Finite')
        X(:, j+1) = exp(-b.*dt).*X(:, j) + fgmcFA(U(:,j), params, Mfft, dt, model, activity);
    end
end

end