function [TCumulants, ECumulants] = computeCumulants(X, params, T, model)
% Cumulants computation in two ways: the first one is the empirical one,
% based on the parameters of the model, while the second one is the
% theoretical one, based on the distribution of the increments of Xt.
% For the theoretical one, first they are computed for the Levy process and
% then used for the OU-Levy or Levy-OU one.
%
% INPUT:
% X:                    discretization of the first input, passed as a vector
% params:               vector of parameters for the specified model
% T:                    maturity
% model:                model selected
%
% OUTPUT:
% TCumulants:           theoretical cumulants for a OU-Levy process
% ECumulants:           empirical cumulants for the simulated OU-Levy process

% Declare for memory allocation.
ECumulants = zeros(1, 4);
TCumulants = zeros(1, 4);

% EMPIRICAL CUMULANTS
ECumulants(1) = mean(X); 
ECumulants(2) = moment(X,2);
ECumulants(3) = moment(X,3);
ECumulants(4) = moment(X,4) - 3*moment(X,2)^2;
ECumulants = ECumulants * 1000;

% THEORETICAL CUMULANTS
switch model
    case {'OU-NTS', 'NTS-OU'}
        % Assign parameters.
        alpha = params(1); b = params(2); sigma = params(3); 
        k = params(4); theta = params(5);
        
        if theta ~= 0 % Asymmetric case of the NTS.
            
            for j = 1:numel(TCumulants)
                if (j == 1)
                    TCumulants(j) = theta;
                elseif (j >= 2)
                    for n = 0:floor(j / 2)
                        TCumulants(j) = TCumulants(j) + factorial(j) / (factorial(n) * factorial(j - 2 * n)) * ...
                            gamma(j - alpha - n) / gamma(1 - alpha) * (k / (1 - alpha))^(j - 1 - n) * theta^(j - 2 * n) * (sigma^2 / 2)^n;
                    end
                end
                if strcmp(model, 'OU-NTS')
                    TCumulants(j) = (1 - exp(-j * b * T)) / (j * b) * TCumulants(j) * 1000;
                elseif strcmp(model, 'NTS-OU')
                    TCumulants(j) = (1 - exp(-j * b * T)) * TCumulants(j) * 1000;
                end
            end 
        else % Symmetric case of the NTS
            for j = 1:numel(TCumulants)
                if mod(j, 2) == 0
                    TCumulants(j) = factorial(j) / factorial(j / 2) * gamma(j / 2 - alpha) / ...
                        gamma(1 - alpha) * (k / (1 - alpha))^(j / 2 - 1) * (sigma^2 / 2)^(j / 2);
                else
                    TCumulants(j) = 0;
                end
                if strcmp(model, 'OU-NTS')
                    TCumulants(j) = (1 - exp(-j * b * T)) / (j * b) * TCumulants(j) * 1000;
                elseif strcmp(model, 'NTS-OU')
                    TCumulants(j) = (1 - exp(-j * b * T)) * TCumulants(j) * 1000;
                end
            end
        end

    case {'OU-TS', 'TS-OU'}
        % Assign parameters.
        alpha = params(1); b = params(2); beta_p = params(3); beta_n = params(4);
        c_p = params(5); c_n = params(6); gamma_c = params(7);
        
        for j = 1:numel(TCumulants)
            if (j == 1)
                TCumulants(j) = gamma_c;
            elseif (j >= 2)
                TCumulants(j) = c_p * beta_p^(alpha - j) * gamma(j - alpha) + ...
                                  (-1)^j*c_n * beta_n^(alpha - j) * gamma(j - alpha);
            end
            if strcmp(model, 'OU-TS')
                TCumulants(j) = (1 - exp(-j * b * T)) / (j * b) * TCumulants(j) * 1000;
            elseif strcmp(model, 'TS-OU')
                TCumulants(j) = (1 - exp(-j * b * T)) * TCumulants(j) * 1000;
            end
        end
end
end