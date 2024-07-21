function [du, a] = extraParamsComputation(model, activity, params, N, dt)
% Proceed to assign parameters, compute analyticity strip (p_n, p_p) and
% compute the optimal step du for the FFT code.
% For the finite activity processes du was found by checking the one that
% gave empirical cumulants closer to the theoretical ones.
%
% Finite activity processes have power decay, infinite activity have
% exponential decay. Special cases are commented.
%
% INPUT:
% model:                model selected
% activity:             model activity, needed for the FA case
% params:               vector of paramteres for the specified model
% N:                    FFT parameter, 2^M
% dt:                   time step
%
% OUTPUT:
% du:                   discretization for the FFT algorithm
% a:                    shift for the integral for CDF reconstruction

alpha = params(1); b = params(2);

switch model
    case {'OU-NTS', 'NTS-OU'}
        sigma = params(3);
        k = params(4); theta = params(5);
    
        A_as = sqrt(theta^2 + (2 * sigma^2 * (1 - alpha)) / k);
        p_n = (theta - A_as) / (sigma^2);
        p_p = (theta + A_as) / (sigma^2);
        
        switch activity
            case 'Finite'
                a = abs(0.25 * max(-p_n, p_p));
                du = 0.0971; % Best for NTS
            case 'Infinite'
                if alpha == 0 % OU-VG case (Variance Gamma, NTS), power decay.
                    a = 0.25 * max(-p_n, p_p);
                    du = 0.0971;
                else
                    a = 0.5 * max(-p_n, p_p);
                    omega = 2 * alpha;
                    if strcmp(model, 'OU-NTS') % OU-NTS GENERAL CASE
                        l = 1 / alpha * ((1 - alpha) / k)^(1 - alpha) * (sigma^2 / 2)^alpha * (1 - exp(-2 * alpha * b * dt)) / (2 * alpha * b);
                    else % NTS-OU GENERAL CASE
                        l = 1 / alpha * ((1 - alpha) / k)^(1 - alpha) * (sigma^2 / 2)^alpha * (1 - exp(-2 * alpha * b * dt));
                    end
                    du = (2 * pi * abs(a) / (l * N^omega))^(1 / (omega + 1));
                end
        end
        
    case {'OU-TS', 'TS-OU'}
        beta_p = params(3); beta_n = params(4); c_p = params(5); c_n = params(6);
    
        p_n = -beta_p;
        p_p = beta_n;

        switch activity
            case 'Finite'
                a = abs(0.25 * max(-p_n, p_p));
                du = 0.1; % Best for TS
            case 'Infinite'
                if alpha == 0 % OU-Gamma case (TS), power decay.
                    a = 0.25 * max(-p_n, p_p);
                    du = 0.0971;
                else
                    a = 0.5 * max(-p_n, p_p);
                    if strcmp(model, 'OU-TS')
                        if alpha == 1 % Other special case for OU-TS, exp decay.
                            omega = 1;
                            l = (c_p + c_n) * pi / 2 * (1 - exp(-b * dt)) / b;
                        else % GENERAL OU-TS CASE.
                            omega = alpha;
                            scale = 1;
                            l = -(c_p + c_n) * gamma(-alpha) * cos(alpha * pi / 2) * (1 - exp(-alpha * b * dt)) / (alpha * b) * scale^alpha;
                        end
                    else % TS-OU
                        if alpha == 1 % Other special case for TS-OU, exp decay.
                            omega = 1;
                            l = (c_p + c_n) * pi / 2 * (1 - exp(-b * dt));
                        else % GENERAL TS-OU CASE.
                            omega = alpha;
                            scale = 1;
                            l = -(c_p + c_n) * gamma(-alpha) * cos(alpha * pi / 2) * (1 - exp(-alpha * b * dt)) * scale^alpha;
                        end
                    end
                    du = (2 * pi * abs(a) / (l * N^omega))^(1 / (omega + 1));
                end
        end
end

end