function Xs = spotSimulationOU(params, N, M, T, seed)
% Simulation of the Spot price with a mean reverting model, using also AV
%
% INPUT:
% mu:                   long mean term
% sigma:                volatility term
% k:                    mean reversion term
% N:                    number of paths
% M:                    number of timesteps
% T:                    last day of Spot trading
% seed:                 (optional) fix the randomness to a specific seed
%
% OUTPUT:
% X:                    logprice from t = 0 to t = T+1
% XAV:                  logprice from t = 0 to t = T+1, AV method

% Set given seed
rng(seed)

dt = T/M; % daily simulation
sigma = params(1); k = params(2); mu = params(3);

% simulate the grid of the exp Xt, the first row is X_0, the last one is X_(T+1) 
X = zeros(N, T+1); 
XAV = zeros(N, T+1);
Z = randn(N, T);

for i=1:T
    X(:, i+1) = X(:,i) + k*(mu-X(:, i) - sigma.^2/(2*k))*dt + sigma*sqrt(dt)*Z(:, i);
    XAV(:, i+1) = XAV(:,i) + k*(mu-XAV(:, i) - sigma.^2/(2*k))*dt - sigma*sqrt(dt)*Z(:, i);
end

Xs = [X; XAV];

end