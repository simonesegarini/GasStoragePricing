function increments = exDeIA(params, dt, nSim, model)
% Compute the increments for TS-OU and OU-TS processes following Algorithm
% 1 and Algorithm 2 of Sabino & Cufaro Petroni 2022.
%
% INPUT:
% params:               vector with the parameters of the model
% dt:                   time step
% nSim:                 number of simulations
% model:                char for model selection
% seed:                 set the seed of the simulation
%
% OUTPUT:
% increments:           increments using exact decomposition method

% Assign parameters.
alpha = params(1); b = params(2); beta_p = params(3);
beta_n = params(4); c_p = params(5); c_n = params(6);
gamma_c = params(7);

% As exDe Algortihm of Sabino & Cufaro Petroni, we decompose 
% X(t) = aX0 + X1 + X2.
% aX0 is already computed in the update rule of exDe.

a = exp(-b*dt);

switch model
    case 'OU-TS'
        % Compute parameters for the poisson rv.
        lamb_p = c_p*beta_p^alpha*gamma(1-alpha) * (1 - a^alpha + a^alpha * log(a^alpha)) / (b*alpha^2 * a^alpha);
        lamb_n = c_n*beta_n^alpha*gamma(1-alpha) * (1 - a^alpha + a^alpha * log(a^alpha)) / (b*alpha^2 * a^alpha);
        
        % Compute X1 by doing a TS simulations.
        X1p = simulationTS_SSR(alpha, beta_p/a, c_p*(1-a^alpha)/(alpha*b), nSim);
        X1n = simulationTS_SSR(alpha, beta_n/a, c_n*(1-a^alpha)/(alpha*b), nSim);

        % Cumulants to compensate.
        ctsCumulants = computeCumulants(0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], dt, 'OU-CTS')/1000;

    case 'TS-OU'
        % Compute parameters for the poisson rv.
        lamb_p = c_p*gamma(1-alpha)*beta_p^alpha/alpha * (1 - a^alpha);
        lamb_n = c_n*gamma(1-alpha)*beta_n^alpha/alpha * (1 - a^alpha);

        % Compute X1 by doing a TS simulations.
        X1p = simulationTS_SSR(alpha, beta_p, c_p*(1-a^alpha), nSim);
        X1n = simulationTS_SSR(alpha, beta_n, c_n*(1-a^alpha), nSim);

        % Cumulants to compensate.
        ctsCumulants = computeCumulants(0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], dt, 'CTS-OU')/1000;
end

% Compute X2 as a compound poisson processes.
X2p = zeros(nSim, 1);
X2n = zeros(nSim, 1);

% Simulate the poisson for the compound poisson.
positivePoisson = poissrnd(lamb_p, nSim, 1);
negativePoisson = poissrnd(lamb_n, nSim, 1);

% Compute max number of jumps, helps to deal with the size of vector Vp Vn.
nP = max(positivePoisson); 
nN = max(negativePoisson);

% Generate n iid rv's with Eq 38 & 42 of Sabino & Cufaro Petroni.
Vp = zeros(nSim, nP);
Vn = zeros(nSim, nN);

switch model
    case 'OU-TS'

        for itP = 1:nP
            Vp(:, itP) = simulationV(params, dt, nSim);
        end
        
        for itN = 1:nN
            Vn(:, itN) = simulationV(params, dt, nSim);
        end

    case 'TS-OU'

        Usp = rand(nSim, nP);
        Usn = rand(nSim, nN);

        Vp = (1 + (a^(-alpha) - 1)./(alpha) .* Usp).^(1/alpha);
        Vn = (1 + (a^(-alpha) - 1)./(alpha) .* Usn).^(1/alpha);
end

% Adjust the values of beta_p and beta_n.
beta_p_hat = beta_p*Vp;
beta_n_hat = beta_n*Vn;

% Generate n independent generalized gamma rvs.
Jip = gamrnd(1-alpha, 1./beta_p_hat);
Jin = gamrnd(1-alpha, 1./beta_n_hat);

% Compute the values of the positive and negative compound poisson for the
% jumps of the process.
for simIt = 1:nSim

    % Positive jumps.
    pPos = positivePoisson(simIt);
    if pPos > 0
        X2p(simIt) = sum(Jip(simIt, 1:pPos));
    end

    % Negative jumps.
    pNeg = negativePoisson(simIt);
    if pNeg > 0
        X2n(simIt) = sum(Jin(simIt, 1:pNeg));
    end
end

% Compute the total increment as return value.
increments = gamma_c*dt + X1p + X2p - X1n - X2n - ctsCumulants(1);
end