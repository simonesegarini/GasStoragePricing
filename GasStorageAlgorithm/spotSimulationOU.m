function Xs = spotSimulationOU(params, N, M, T, seed)
% Simulation of the Spot price with a mean reverting model, using also AV
%
% INPUT:
% params:               parameters for the OU algorithm
% N:                    number of paths
% M:                    number of timesteps
% T:                    last day of Spot trading
% seed:                 (optional) fix the randomness to a specific seed
%
% OUTPUT:
% Xs:                   logprice from t = 0 to t = T+1

% Set given seed.
rng(seed)

dt = T/M; % daily simulation
sigma = params(1); k = params(2); A = params(3); B = params(4);
mu = @(t) A.*cos(2.*pi.*t./T) + B;

% Simulate the grid of the exp Xt, the first row is X_0, the last one is
% X_(T+1).
X = zeros(N, M+1); 
XAV = zeros(N, M+1);
Z = randn(N, M);

% Forward procedure using the discreziation for the SDE.
for j=1:M
    X(:, j+1) = X(:,j) + k*(mu(j)-X(:, j) - sigma.^2/(2*k))*dt + sigma*sqrt(dt)*Z(:, j);
    XAV(:, j+1) = XAV(:,j) + k*(mu(j)-XAV(:, j) - sigma.^2/(2*k))*dt - sigma*sqrt(dt)*Z(:, j);
end

Xs = [X; XAV];

end