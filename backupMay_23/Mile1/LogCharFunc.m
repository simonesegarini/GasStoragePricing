function values = LogCharFunc(us, t, params, model, activity)
% Log Charachteristic Functionof a Levy process, return the values given a
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

% TS DO THE BILATERAL CASE
% DECLARE HANDLES OUTSIDE OF FOR LOOP

% Preallocate for speed
values = zeros(size(us));

switch model
    case 'OU-NTS'   % DONE VVV
        % assign parameters
        alpha = params(1); b = params(2); sigma = params(3);
        k = params(4); theta = params(5);
        if strcmp(activity, 'Infinite')
            % CHECK THIS CASE
            if alpha == 0 % VG CASE, POWER LAW FOR AS BEHAVIOUR, B.15/B.16
                integrand = @(z, u) log(0.5*sigma.^2.*u.^2.*z.^2.*k - 1i.*theta.*k.*u.*z + 1)./z;
    
                for i=1:numel(us)
                    z_values = linspace(exp(-b*t), 1, 1000);
                    integral_value = trapz(z_values, integrand(z_values, us(i)));
        
                    values(i) = -1/(k*b)*integral_value;
                end
            else % GENERAL CASE FOR OU-NTS, EXP LAW FOR ALPHA > 0, POWER ALPHA < 0, .42
                integrand = @(z, u) (0.5.*sigma.^2.*u.^2.*z.^2 - 1i*theta.*u.*z + (1-alpha)./k).^alpha./z;
    
                for i=1:numel(us)
                    z_values = linspace(exp(-b*t), 1, 1000); % Adjust the number of points as needed
                    integral_value = trapz(z_values, integrand(z_values, us(i)));
                    
                    values(i) = (1-alpha)./(k*alpha).*t - (((1-alpha)/k).^(1-alpha))./(alpha*b).*integral_value;
                end
            end
        elseif strcmp(activity, 'Finite')
            aux_values = zeros(size(values));
            integrand = @(z, u) (0.5.*sigma.^2.*u.^2.*z.^2 - 1i*theta.*u.*z + (1-alpha)./k).^alpha./z;
    
            for i=1:numel(us)    %.41 
                z_values = linspace(exp(-b*t), 1, 1000);
                integral_value = trapz(z_values, integrand(z_values, us(i)));
                
                aux_values(i) = ((1-alpha)./k).^(-alpha)./(b*t).*integral_value;
            end
            lamb = (1-alpha)./(k*abs(alpha));
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end
    case 'NTS-OU'
        alpha = params(1); b = params(2); sigma = params(3);
        k = params(4); theta = params(5);
        if strcmp(activity, 'Infinite')
            % GENERAL CASE FOR NTS-OU, EXP LAW FOR ALPHA > 0, POWER ALPHA < 0, .43
            values = (1-alpha)/(k*alpha) .* ((1-(1i.*k)./(1-alpha) .* (theta.*us.*exp(-b*t) + ...
                1i.*sigma^2.*us.^2.*exp(-2*b*t)./2)).^alpha - (1-(1i.*k)./(1-alpha) .* (theta.*us + ...
                1i.*sigma^2.*us.^2./2)).^alpha);

        elseif strcmp(activity, 'Finite')
            aux_values = 1/(2*b*t).*log((0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*us.*exp(b*t) + exp(2*b*t))./...
                    (0.5*sigma^2.*us.^2.*k - 1i.*theta.*k.*u + 1));
            lamb = 2*b/k;
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end

    case'OU-TS' %FIX XXX
        % assign parameters
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);
        if strcmp(activity, 'Infinite')
            if alpha == 0 % GAMMA CASE, POWER LAW FOR AS BEHAVIOUR, B.3/B.4, USE ALG1
                integrand = @(z, u) log(1-1i.*u./z)./z;
                
                for i=1:numel(us)
                    z_values = linspace(beta_p, beta_p*exp(b*t),1000);
                    integral_values = trapz(z_values, integrand(z_values, us(i)));
        
                    values(i) = 1i.*us(i).*(1-exp(-b.*t))./b.*(gamma_c - c_p./beta_p) - ...
                        c_p./b.*integral_values(i);
                end
            elseif alpha == 1 % OTH SPECIAL TS CASE, EXPONENTIAL LAW FOR AS BEHAVIOUR, B.9/B.10
                integrand = @(z, u) (1-1i.*u./z).*log(1-1i.*u./z)./z;
    
                for i=1:numel(us)
                    z_values = linspace(beta_p, beta_p*exp(b*t),1000);
                    integral_values = trapz(z_values, integrand(z_values, us(i)));
        
                    values(i) = 1i.*us(i).*(1-exp(-b.*t))./b.*(gamma_c + c_p) + ...
                        c_p.*beta_p./b.*integral_values;
                end
            else % GENERAL CASE FOR TS, EXP LAW FOR AS BEHAVIOUR .15
                integrand1 = @(z, u) ((z-1i.*u).^alpha)./(z.^(alpha+1));
                integrand2 = @(z, u) ((z+1i.*u).^alpha)./(z.^(alpha+1));
    
                for i=1:numel(us)
                    z_values1 = linspace(beta_p, beta_p.*exp(b.*t),1000);
                    z_values2 = linspace(beta_n, beta_n.*exp(b.*t),1000);
                    integral_value1 = trapz(z_values1, integrand1(z_values1, us(i)));
                    integral_value2 = trapz(z_values2, integrand2(z_values2, us(i)));
        
                    values(i) = 1i.*us(i).*(1-exp(b.*t))./b.*gamma_c + ...
                        c_p.*beta_p.^alpha.*gamma(-alpha)./b.*(integral_value1 ...
                        - b.*t + alpha./beta_p.*1i.*us(i).*(1-exp(-b.*t))) + ...
                        c_n.*beta_n.^alpha.*gamma(-alpha)./b.*(integral_value2 ...
                        - b.*t - alpha./beta_n.*1i.*us(i).*(1-exp(-b.*t)));
                end
            end
        elseif strcmp(activity, 'Finite')
            aux_values = zeros(size(values));
            integrand = @(z, u) (z-1i.*u).^alpha./z.^(alpha+1);
    
            for i=1:numel(us)    %BISOGNA IMPLEMENTARE ANCHE LA PARTE DI C_N .38
                z_values = linspace(beta_p, beta_p*exp(b*t));
                integral_value = trapz(z_values, integrand(z_values, us(i)));
                
                aux_values(i) = 1./(b*t).*integral_value;
            end
            lamb = c_p.*beta_p*gamma(-alpha);
            values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
        end
    case 'TS-OU'
        % assign parameters
        alpha = params(1); b = params(2); beta_p = params(3);
        beta_n = params(4); c_p = params(5); c_n = params(6);
        gamma_c = params(7);
        if strcmp(activity, 'Infinite')
            if alpha == 1 % OTH SPECIAL TS CASE, EXPONENTIAL LAW FOR AS BEHAVIOUR, B.9/B.10
                values = 1i.*us.*(1-exp(-b*t)).*(gamma_c+c_p) + ...
                    c_p*beta_p*((1- 1i.*us./beta_p).*log(1- 1i.*us./beta_p) - ...
                    (1-1i.*us.*exp(-b*t)./beta_p).*log(1-1i.*us.*exp(-b*t)./beta_p));
            else % GENERAL CASE FOR TS, EXP LAW FOR AS BEHAVIOUR .16
                values = 1i.*us.*(1-exp(-b*t)).*...
                    (gamma_c - c_p.*gamma(1-alpha)).*beta_p.^(alpha-1) + c_n.*gamma(1-alpha).*beta_n.^(alpha-1) +...
                    c_p.*gamma(-alpha).*((beta_p-1i.*us).^alpha - (beta_p-1i.*us.*exp(-b*t)).^alpha) + ...
                    c_n.*gamma(-alpha).*((beta_n-1i.*us).^alpha - (beta_n-1i.*us.*exp(-b*t)).^alpha);
            end
        elseif strcmp(activity, 'Finite')
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