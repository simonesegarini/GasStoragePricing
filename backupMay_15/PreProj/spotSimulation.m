function [X, XAV] = spotSimulation(mu, sigma, k, M, T, seed)
% Simulation of the Spot price with a mean reverting model, using also AV
%
% INPUT:
% mu:                   long mean term
% sigma:                volatility term
% k:                    mean reversion term
% M:                    number of paths
% T:                    last day of Spot trading
% seed:                 (optional) fix the randomness to a specific seed
%
% OUTPUT:
% X:                    logprice from t = 0 to t = T+1
% XAV:                  logprice from t = 0 to t = T+1, AV method

% check if the user gives a seed
if nargin == 6
    rng(seed)
end

dt = 1; % daily simulation

% simulate the grid of the exp Xt, the first row is X_0, the last one is X_(T+1) 
X = zeros(M, T+1); 
XAV = zeros(M, T+1);
Z = randn(M, T);

for i=1:T
    X(:, i+1) = X(:,i) + k*(mu(i)-X(:, i) - sigma.^2/(2*k))*dt + sigma*sqrt(dt)*Z(:, i);
    XAV(:, i+1) = XAV(:,i) + k*(mu(i)-XAV(:, i) - sigma.^2/(2*k))*dt - sigma*sqrt(dt)*Z(:, i);
end

end