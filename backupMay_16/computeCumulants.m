function [TCumulants, ECumulants] = computeCumulants(X, params, T, model)
% Cumulants computation in two ways: the first one is the theoretic one,
% based on the parameters of the model, while teh second one is the
% empirical one, based on the distribution of the increments of Xt
%
% INPUT:
% Xt:                   discretization of the first input, passed as a vector
% params:               vector of paramteres for the specified model
% T:                    maturity
% model:                model selected
%
% OUTPUT:
% TCumulants:           theoretical cumulants for a OU-Levy process
% ECumulants:           empirical cumulants for the simulated OU-Levy process

% Declare for memory allocation
ECumulants = zeros(1,4);
TCumulants = zeros(1,4);

% EMPIRICAL CUMULANTS
ECumulants(1) = mean(X); 
ECumulants(2) = var(X);
ECumulants(3) = skewness(X) * sqrt(ECumulants(2))^3 - 3*(ECumulants(2) + ...
    ECumulants(1)^2)*ECumulants(1) + 2*ECumulants(1)^3; 
ECumulants(4) = kurtosis(X)*ECumulants(2)^2 - 4*skewness(X)*sqrt(ECumulants(2))^3*ECumulants(1) - ...
    3*(ECumulants(2)+ECumulants(1)^2)^2 + 12*(ECumulants(2) + ECumulants(1)^2)*ECumulants(1)^2 - 6*ECumulants(1)^4;
ECumulants = ECumulants*1000;

% THEORETICAL CUMULANTS
if strcmp(model, 'NTS')
    % Assign paramteres
    alpha = params(1); b = params(2); sigma = params(3); 
    k = params(4); theta = params(5);
    
    % EXTRA: can be implementend the formula that simplifies the symm case, .23
    if theta ~= 1
        for kit = 1:numel(TCumulants)
            if (kit == 1)
                TCumulants(kit) = theta;
            elseif (kit >= 2)
                for n = 0:floor(kit/2)
                    TCumulants(kit) = TCumulants(kit) + factorial(kit)/(factorial(n)*factorial(kit-2*n)) * ...
                        gamma(kit-alpha-n)/gamma(1-alpha) * (k/(1-alpha))^(kit-1-n)*theta^(kit-2*n)*(sigma^2/2)^n;
                end
            end
    
            TCumulants(kit) = (1-exp(-kit*b*T))/(kit*b)*TCumulants(kit)*1000;
        end
    elseif theta == 0 % If NTS is symmetric we have the simple expression
        for kit = 1:numel(TCumulants)
            if mod(kit, 2) == 0
                ECumulants(kit) = factorial(kit)/factorial(k/2) * gamma(kit/2-alpha) / ...
                    gamma(1-alpha) * (k/(1-alpha))^(kit/2-1) * (sigma^2/2)^(kit/2);
            else
                ECumulants(kit) = 0;
            end
        end
    end

elseif strcmp(model, 'TS')
    % Assign parameters
    alpha = params(1); b = params(2); beta_p = params(3); beta_n = params(4);
    c_p = params(5); c_n = params(6); gamma_c = params(7);
    
    for kit = 1:numel(TCumulants)
        if (kit == 1)
            TCumulants(kit) = gamma_c;
        elseif (kit >= 2)
            TCumulants(kit) = c_p*beta_p^(alpha-kit)*gamma(kit-alpha) +...
                c_n*beta_n^(alpha-kit)*gamma(kit-alpha);
        end

        TCumulants(kit) = (1-exp(-kit*b*T))/(kit*b)*TCumulants(kit)*1000;
    end
end
end