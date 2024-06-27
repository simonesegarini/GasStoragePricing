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
    case 'OU-NTS'   % DONE VVV
        % assign parameters
        alpha = params(1); b = params(2); sigma = params(3);
        k = params(4); theta = params(5);
        if strcmp(activity, 'Infinite')
            % CHECK THIS CASE
            if alpha == 0 % VG CASE, POWER LAW FOR AS BEHAVIOUR, B.15/B.16
                integral_values = integral(@(z) log(0.5*sigma.^2.*us.^2.*z.^2.*k - 1i.*theta.*k.*us.*z + 1)./z, exp(-b*t), 1, 'ArrayValued', true);
                values = -1/(k*b)*integral_values;

            else % GENERAL CASE FOR OU-NTS, EXP LAW FOR ALPHA > 0, POWER ALPHA < 0, .42
                integral_values = integral(@(z) (0.5.*sigma.^2.*us.^2.*z.^2 - 1i*theta.*us.*z + (1-alpha)./k).^alpha./z, exp(-b*t), 1, 'ArrayValued', true);
                values = (1 - alpha) ./ (k * alpha) .* t - (((1 - alpha) / k).^(1 - alpha)) ./ (alpha * b) .* integral_values;
            end

        elseif strcmp(activity, 'Finite')
            integral_values = integral(@(z) (0.5.*sigma.^2.*us.^2.*z.^2 - 1i*theta.*us.*z + (1-alpha)./k).^alpha./z, exp(-b*t), 1, 'ArrayValued', true);
            aux_values = ((1-alpha)./k).^(-alpha)./(b*t).*integral_values;

            lamb = (1-alpha)./(k*abs(alpha));
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end
    case 'NTS-OU' %DONE VVV
        alpha = params(1); b = params(2); sigma = params(3);
        k = params(4); theta = params(5);
        if strcmp(activity, 'Infinite')
            % GENERAL CASE FOR NTS-OU, EXP LAW FOR ALPHA > 0, POWER ALPHA < 0, .43
            values = psiXNTS(us, params) - psiXNTS(us*exp(-b*t), params);

        elseif strcmp(activity, 'Finite')
            aux_values = 1/(2*b*t).*log((0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*us.*exp(b*t) + exp(2*b*t))./...
                    (0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*u + 1));
            lamb = 2*b/k;
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end

    case'OU-TS' %DONE VVV
        % assign parameters
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);
        if strcmp(activity, 'Infinite')
            if alpha == 0 % GAMMA CASE, POWER LAW FOR AS BEHAVIOUR, B.3/B.4
                integral_values1 = integral(@(z) log(1-1i.*us./z)./z, beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) log(1+1i.*us./z)./z, beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(-b.*t))./b.*(gamma_c - c_p./beta_p + c_n./beta_n) - ...
                                         c_p./b.*integral_values1 - c_n./b.*integral_values2;

            elseif alpha == 1 % OTH SPECIAL TS CASE, EXPONENTIAL LAW FOR AS BEHAVIOUR, B.9/B.10
                integral_values1 = integral(@(z) (1-1i.*us./z).*log(1-1i.*us./z)./z, beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) (1-1i.*us./z).*log(1-1i.*us./z)./z, beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(-b.*t))./b.*(gamma_c + c_p + c_n) + ...
                                        c_p.*beta_p./b.*integral_values1.' + c_n.*beta_n./b.*integral_values2.';

            else % GENERAL CASE FOR TS, EXP LAW FOR AS BEHAVIOUR .15
                integral_values1 = integral(@(z) (z-1i*us).^alpha./z.^(alpha+1), beta_p, beta_p*exp(b*t), 'ArrayValued', true);
                integral_values2 = integral(@(z) (z+1i*us).^alpha./z.^(alpha+1), beta_n, beta_n*exp(b*t), 'ArrayValued', true);

                values = 1i.*us.*(1-exp(b.*t))./b.*gamma_c + ...
                                        c_p.*beta_p.^alpha.*gamma(-alpha)./b.*(integral_values1 ...
                                        - b.*t + alpha./beta_p.*1i.*us.*(1-exp(-b.*t))) + ...
                                        c_n.*beta_n.^alpha.*gamma(-alpha)./b.*(integral_values2 ...
                                        - b.*t - alpha./beta_n.*1i.*us.*(1-exp(-b.*t)));
            end
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
    case 'TS-OU' %%% XXXX FIX
        % assign parameters
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);
        %%%%%%%%%%%%% CN
        if strcmp(activity, 'Infinite')
            if alpha == 1 % OTH SPECIAL TS CASE, EXPONENTIAL LAW FOR AS BEHAVIOUR, B.9/B.10
                %FIX FIX FIX
                values = 1i.*us.*(1-exp(-b*t)).*(gamma_c+c_p) + ...
                    c_p*beta_p*((1- 1i.*us./beta_p).*log(1- 1i.*us./beta_p) - ...
                    (1-1i.*us.*exp(-b*t)./beta_p).*log(1-1i.*us.*exp(-b*t)./beta_p));

            else % GENERAL CASE FOR TS, EXP LAW FOR AS BEHAVIOUR .16
                % DONE DONE DONE
                values = psiXTS(us, params) - psiXTS(us*exp(-b*t), params);
            end

        elseif strcmp(activity, 'Finite') %FIX FIX FIX
            % GAMMA CASE B.6
            aux_values = zeros(size(values));
            integrand = @(z, u) 1/(z-1i.*u);
    
            for i=1:numel(us)    %BISOGNA IMPLEMENTARE ANCHE LA PARTE DI C_N .38
                z_values = linspace(beta_p, beta_p*exp(b*t));
                integral_value = trapz(z_values, integrand(z_values, us(i)));
                
                aux_values(i) = 1./(b*t).*integral_value;
            end
            lamb = c_p.*b;
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end
end
end