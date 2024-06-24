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

% Assign parameters in order to compute lambda to apply Algorithm 2 for finite
% activity processes.
switch model
    case 'OU-NTS'
        alpha = params(1);k = params(4);
    
        lamb = (1-alpha)./(k*abs(alpha));
    case 'NTS-OU'
        % We have the special case of VG-OU that is FA but the parameters
        % lambda is the same, the part that change is in LogCharFunc
        % function.
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

        lamb_p = c_p*b; 
        lamb_n = c_n*b;
        lamb = lamb_p + lamb_n;
end

% Check where we have to compute the increment, i.e. where we have jumps in
% the finite activity model.
Bt = rand(size(U)) < 1-exp(-lamb*dt);
idxs = find(Bt == 1);

% Compute increments by calling the function for the inifinite activity
% processes just where we have a jump.
increments = zeros(size(U));
increments(idxs) = fgmcIA(U(idxs), params, M, dt, model, activity);
end