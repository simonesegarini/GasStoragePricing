function se = standardError(x)
% Compute the standrad error of the given vector, used as a metric to check
% that increasing the number of simulations works.
%
% INPUT:
% x:                input vector of prices
%
% OUTPUT:
% se:               stadard error

% Compute the standard deviation of the vector x
sigma = std(x);

% Compute the number of observations in the vector x
n = length(x);

% Compute the standard error
se = sigma / sqrt(n);
end
