function increments = exDeFA(params, dt, nSim, model)
% Compute the increments for OU-TS following Sabino 2022 and for TS-OU
% follwing Sabino & Cufaro Petroni 2021.
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

% As exDe Algortihm of Sabino, we decompose X(t) = aX0 + X1.
% aX0 is already computed in the update rule of exDe.

switch model
    case 'OU-TS'
        % Compute parameters for the poisson rv.
        lamb_p = c_p*gamma(-alpha)*beta_p^alpha;
        lamb_n = c_n*gamma(-alpha)*beta_n^alpha;

        % Cumulants to compensate.
        ctsCumulants = computeCumulants(0, [alpha, b, beta_p, beta_n, c_p, c_n, gamma_c], dt, 'OU-CTS')/1000;

        % Compute X1 as a compound poisson processes.
        X1p = zeros(nSim, 1);
        X1n = zeros(nSim, 1);
        
        % Simulate the poisson for the compound poisson.
        positivePoisson = poissrnd(lamb_p*dt, nSim, 1);
        negativePoisson = poissrnd(lamb_n*dt, nSim, 1);
        
        % Compute max number of jumps, helps to deal with the size of vector Usp Usn.
        nP = max(positivePoisson); 
        nN = max(negativePoisson);
        
        % Generate n iid uniform rvs.
        Usp = rand(nSim, nP);
        Usn = rand(nSim, nN);
        
        beta_p_hat = beta_p.*exp(b*Usp*dt);
        beta_n_hat = beta_n.*exp(b*Usn*dt);
        
        % Generate n independent generalized gamma rvs.
        Jip = gamrnd(-alpha, 1./beta_p_hat);
        Jin = gamrnd(-alpha, 1./beta_n_hat);
        
        % Compute the values of the positive and negative compound poisson for the
        % jumps of the process.
        for simIt = 1:nSim
        
            % Positive jumps.
            pPos = positivePoisson(simIt);
            if pPos > 0
                X1p(simIt) = sum(Jip(simIt, 1:pPos));
            end
        
            % Negative jumps.
            pNeg = negativePoisson(simIt);
            if pNeg > 0
                X1n(simIt) = sum(Jin(simIt, 1:pNeg));
            end
        end

    % The TS-OU process is a finite activity only if alpha = 0 i.e. it is a gamma-OU.
    case 'TS-OU' 
        % Sabino scrive gamma-ou(k, lambda, beta), k è b

        a = exp(-b*dt);
        alphap = lamb_p/k;
        alphan = lamb_n/k;

%         % Generate Polya random variables.
%         Bp = simulationPolya(alphap, 1-a, nSim);
%         Bn = simulationPolya(alphan, 1-a, nSim);
% 
%         % Generate the increments as Erlang random variables.
%         X1p = simulationErlang(Bp, betap/a, nSim);
%         X1n = simulationErlang(Bn, betan/a, nSim);

        
        % shape parameter == alpha;
        % rate = beta;
        Bp = nbinrnd(alphap, a,  size(nSim));
        Bn = nbinrnd(alphan, a,  size(nSim));
        
        X1p = gamrnd(Bp, a/betap, size(nSim));
        X1n = gamrnd(Bn, a/betan, size(nSim));
end

% Compute the total increment as return value.
% increments = gamma_c*dt + X1p - X1n - cumulants_p(1) + cumulants_n(1);
increments = gamma_c*dt + X1p - X1n - ctsCumulants(1);
end
