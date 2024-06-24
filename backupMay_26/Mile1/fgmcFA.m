function increments = fgmcFA(U, params, M, dt, model, activity)
% FGMC method for Finite Activity processes
%
% INPUT:
% U:                    matrix with the simulated Uniform RV
% params:               vector of paramteres for the specified model
% M:                    parameter for FFT
% dt:                   time step
% model:                model selected
% activity:             model activity, needed for the FA case
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy/Levy-OU


% TS DO THE BILATERAL CASE
% Assign parameters and compute analyticity strip
switch model
    case 'OU-NTS'
        alpha = params(1);k = params(4);
    
        lamb = (1-alpha)./(k*abs(alpha));
    case 'NTS-OU'
        b = params(2); k = params(4);

        lamb = 2*b/k;
    case 'OU-TS'
        alpha = params(1); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
    
        lamb_p = c_p*beta_p^alpha*gamma(-alpha);
        lamb_n = c_n*beta_n^alpha*gamma(-alpha);
        lamb = lamb_p+lamb_n;
    case 'TS-OU'
        b = params(2); c_p = params(5); c_n = params(6);

        lamb = cp*b;
end

Bt = rand(size(U)) < 1-exp(-lamb*dt);
idxs = find(Bt == 1);

increments = zeros(size(U));
increments(idxs) = fgmcIA(U(idxs), params, M, dt, model, activity);
end