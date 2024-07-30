function values = LogCharFunc(us, t, params, model, activity)
% Log Charachteristic Function of a Levy process, return the values given a
% vector of u and a specified time t
%
% INPUT:
% u:                    discretization of the first input, passed as a vector
% t:                    time step
% params:               vector of paramteres for the specified model
% model:                model selected
% activity:             model activity, needed for the FA case
%
% OUTPUT:
% values:               LCF computed in u and in t

% Preallocate for speed.
values = zeros(size(us));

switch model
    case 'OU-NTS'
        % Assign parameters.
        alpha = params(1); b = params(2); sigma = params(3);
        k = params(4); theta = params(5);

        % Infinite activity case.
        if strcmp(activity, 'Infinite')

            if alpha == 0 % OU-VG case.
                integral_values = integral(@(z) log(0.5*sigma.^2.*us.^2.*z.^2.*k - 1i.*theta.*k.*us.*z + 1)./z, exp(-b*t), 1, 'ArrayValued', true);
                values = -1/(k*b)*integral_values;

            else % General OU-NTS case.
                integral_values = integral(@(z) (0.5.*sigma.^2.*us.^2.*z.^2 - 1i*theta.*us.*z + (1-alpha)./k).^alpha./z, exp(-b*t), 1, 'ArrayValued', true);
                values = (1 - alpha) ./ (k * alpha) .* t - (((1 - alpha) / k).^(1 - alpha)) ./ (alpha * b) .* integral_values;
            end

        % Finite activity case, use Algorithm 2.
        elseif strcmp(activity, 'Finite')
            integral_values = integral(@(z) (0.5.*sigma.^2.*us.^2.*z.^2 - 1i*theta.*us.*z + (1-alpha)./k).^alpha./z, exp(-b*t), 1, 'ArrayValued', true);
            aux_values = ((1-alpha)./k).^(-alpha)./(b*t).*integral_values;

            lamb = (1-alpha)./(k*abs(alpha));
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end

    case 'NTS-OU'
        % Assign parameters.
        b = params(2); sigma = params(3);
        k = params(4); theta = params(5);

        % Infinite activity case, also the general one.
        if strcmp(activity, 'Infinite')
            values = psiXNTS(us, params) - psiXNTS(us*exp(-b*t), params);

        % VG-OU case. Here is implemented with the direct formula, can be
        % computed also with the decomposed formula by using psiXNTS.
        elseif strcmp(activity, 'Finite')
            aux_values = 1/(2*b*t).*log((0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*us.*exp(b*t) + exp(2*b*t))./...
                    (0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*us + 1));
            lamb = 2*b/k;
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end

    case'OU-TS' 
        % Assign parameters.
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);

        % Infinite activity case.
        if strcmp(activity, 'Infinite')
 
            if alpha == 0 % OU-Gamma case.
                integral_values1 = integral(@(z) log(1-1i.*us./z)./z, beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) log(1+1i.*us./z)./z, beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(-b.*t))./b.*(gamma_c - c_p./beta_p + c_n./beta_n) - ...
                                         c_p./b.*integral_values1 - c_n./b.*integral_values2;

            elseif alpha == 1 % Second OU-TS special case.
                integral_values1 = integral(@(z) (1-1i.*us./z).*log(1-1i.*us./z)./z, beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) (1-1i.*us./z).*log(1-1i.*us./z)./z, beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(-b.*t))./b.*(gamma_c + c_p + c_n) + ...
                                        c_p.*beta_p./b.*integral_values1.' + c_n.*beta_n./b.*integral_values2.';
            
            else % General OU-TS case.
                integral_values1 = integral(@(z) (z-1i*us).^alpha./z.^(alpha+1), beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) (z+1i*us).^alpha./z.^(alpha+1), beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(b.*t))./b.*gamma_c + ...
                                        c_p.*beta_p.^alpha.*gamma(-alpha)./b.*(integral_values1 ...
                                        - b.*t + alpha./beta_p.*1i.*us.*(1-exp(-b.*t))) + ...
                                        c_n.*beta_n.^alpha.*gamma(-alpha)./b.*(integral_values2 ...
                                        - b.*t - alpha./beta_n.*1i.*us.*(1-exp(-b.*t)));
            end

        % Finite activity case, use Algorithm 2. 
        % (Try to figure out how to use the mu(t) parameter)
        elseif strcmp(activity, 'Finite')
            integral_values_p = integral(@(z) (z-1i.*us).^alpha./(z.^(alpha+1)), beta_p, beta_p*exp(b*t), 'ArrayValued', true);
            integral_values_n = integral(@(z) (z+1i.*us).^alpha./(z.^(alpha+1)), beta_n, beta_n*exp(b*t), 'ArrayValued', true);

            phiJ_p = 1./(b*t).*integral_values_p;
            phiJ_n = 1./(b*t).*integral_values_n;

            lamb_p = c_p*beta_p^alpha*gamma(-alpha);
            lamb_n = c_n*beta_n^alpha*gamma(-alpha);
            lamb = lamb_p+lamb_n;
            phiJ = lamb_p/lamb.*phiJ_p  + lamb_n/lamb.*phiJ_n; 
            values = log((exp(lamb*t.*phiJ) - 1)./(exp(lamb*t) - 1));
        end

    case 'TS-OU' 
        % Assign parameters.
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);

        % Infinite activity case.
        if strcmp(activity, 'Infinite')

            if alpha == 1 % Second TS-OU special case. (Need to be fixed cause c_n is missing)
                values = 1i.*us.*(1-exp(-b*t)).*(gamma_c+c_p) + ...
                    c_p*beta_p*((1- 1i.*us./beta_p).*log(1- 1i.*us./beta_p) - ...
                    (1-1i.*us.*exp(-b*t)./beta_p).*log(1-1i.*us.*exp(-b*t)./beta_p));

            else % General case. (This one correctly implemented)
                values = psiXTS(us, params) - psiXTS(us*exp(-b*t), params);
            end
        
        % Gamma-OU case. (Need to be fixed cause c_n is missing)
        elseif strcmp(activity, 'Finite') 
            integral_values = integral(@(z) 1./(z-1i.*us), beta_p, beta_p*exp(b*t), 'ArrayValued', true);
            
            aux_values = 1./(b*t).*integral_values;
            lamb = c_p.*b;
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end
     
end
end