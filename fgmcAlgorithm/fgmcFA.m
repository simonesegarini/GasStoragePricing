function increments = fgmcFA(xgrid_hat, CDF_hat, U, dt, params, model)
% FGMC method for Finite Activity processes.
%
% INPUT:
% xgrid_hat:            xgrid of reconstructed cdf
% CDF_hat               ygrid of reconstructed cdf
% U:                    matrix with the simulated Uniform RV
% dt:                   time step
% params:               vector of paramteres for the specified model
% model:                model selected
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy/Levy-OU

% Assign parameters in order to compute lambda to apply Algorithm 2 for finite
% activity processes.
switch model
    case 'OU-NTS'
        alpha = params(1);k = params(4);
    
        lamb = (1-alpha)./(k*abs(alpha));
        mu = 0;
    case 'NTS-OU'
        b = params(2); k = params(4);

        lamb = 2*b/k;
        mu = 0;
    case 'OU-TS'
        alpha = params(1); b = params(2); beta_p = params(3); 
        beta_n = params(4); c_p = params(5); c_n = params(6); gamma_c = params(7);
    
        lamb_p = c_p*beta_p^alpha*gamma(-alpha);
        lamb_n = c_n*beta_n^alpha*gamma(-alpha);
        lamb = lamb_p+lamb_n;
        mu = (1-exp(-b.*dt))./b.*(gamma_c + lamb_p.*alpha./beta_p - lamb_n.*alpha./beta_n);

    case 'TS-OU'
        b = params(2); beta_p = params(3); 
        beta_n = params(4); c_p = params(5); c_n = params(6); gamma_c = params(7);

        lamb_p = c_p*b; 
        lamb_n = c_n*b;
        lamb = lamb_p + lamb_n;
        mu = (1-exp(-b*dt)) * (gamma_c - c_p/beta_p + c_n/beta_n);
end

% Check where we have to compute the increment, i.e. where we have jumps in
% the finite activity model.
Bt = rand(size(U)) < 1-exp(-lamb*dt);
idxs = find(Bt == 1);

% Compute increments by calling the function for the inifinite activity
% processes just where we have a jump.
increments = zeros(size(U));
increments(idxs) = fgmcIA(xgrid_hat, CDF_hat, U(idxs));
increments = increments + mu;
end