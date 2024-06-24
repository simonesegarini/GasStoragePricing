function values = LogCharFunc(u, t, params, model, activity)
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

% Preallocate for speed
values = zeros(size(u));
if nargin == 4
    activity = 'Infinite';
end

if strcmp(model, 'OU-NTS')
    % assign parameters
    alpha = params(1); b = params(2); sigma = params(3);
    k = params(4); theta = params(5);
    if strcmp(activity, 'Infinite')
        if alpha == 0 % VG CASE, POWER LAW FOR AS BEHAVIOUR, B.15/B.16
            for i=1:numel(u)
                integrand = @(z) log(0.5*sigma.^2.*u(i).^2.*z.^2.*k - 1i.*theta.*k.*u(i).*z + 1)./z;
                z_values = linspace(exp(-b*t), 1, 1000);
                integral_value = trapz(z_values, integrand(z_values));
    
                values(i) = -1/(k*b)*integral_value;
            end
        else % GENERAL CASE FOR NTS, EXP LAW FOR ALPHA > 0, POWER ALPHA < 0, .42
            for i=1:numel(u)
                integrand = @(z) (0.5.*sigma.^2.*u(i).^2.*z.^2 - 1i*theta.*u(i).*z + (1-alpha)./k).^alpha./z;
                z_values = linspace(exp(-b*t), 1, 1000); % Adjust the number of points as needed
                integral_value = trapz(z_values, integrand(z_values));
                
                values(i) = (1-alpha)./(k*alpha).*t - (((1-alpha)/k).^(1-alpha))./(alpha*b).*integral_value;
            end
        end
    elseif strcmp(activity, 'Finite')
        aux_values = zeros(size(values));
        for i=1:numel(u)    %.41
            integrand = @(z) (0.5.*sigma.^2.*u(i).^2.*z.^2 - 1i*theta.*u(i).*z + (1-alpha)./k).^alpha./z;
            z_values = linspace(exp(-b*t), 1, 1000);
            integral_value = trapz(z_values, integrand(z_values));
            
            aux_values(i) = ((1-alpha)./k).^(-alpha)./(b*t).*integral_value;
        end
        lamb = (1-alpha)./(k*abs(alpha));
        values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
    end

elseif strcmp(model, 'OU-TS')
    % assign parameters
    alpha = params(1); b = params(2); beta_p = params(3);
    beta_n = params(4); c_p = params(5); c_n = params(6);
    gamma_c = params(7);
    if strcmp(activity, 'Infinite')
        if alpha == 0 % GAMMA CASE, POWER LAW FOR AS BEHAVIOUR, B.3/B.4, USE ALG1
            for i=1:numel(u)
                integrand = @(z) log(1-1i.*u(i)/z)./z;
                z_values = linspace(beta_p, beta_p*exp(b*t),1000);
                integral_values = trapz(z_values, integrand(z_values));
    
                values(i) = 1i.*u.*(1-exp(-b.*t))./b.*(gamma_c - c_p./beta_p) - ...
                    c_p./b.*integral_values(i);
            end
        elseif alpha == 1 % OTH SPECIAL TS CASE, EXPONENTIAL LAW FOR AS BEHAVIOUR, B.9/B.10
            for i=1:numel(u)
                integrand = @(z) (1-1i.*u(i)./z).*log(1-1i.*u(i)/z)./z;
                z_values = linspace(beta_p, beta_p*exp(b*t),1000);
                integral_values = trapz(z_values, integrand(z_values));
    
                values(i) = 1i.*u(i).*(1-exp(-b.*t))./b.*(gamma_c + c_p) + ...
                    c_p.*beta_p./b.*integral_values;
            end
        else % GENERAL CASE FOR TS, EXP LAW FOR AS BEHAVIOUR .38
            for i=1:numel(u)
                integrand1 = @(z) (z-1i.*u(i)).^alpha./z.^(alpha+1);
                integrand2 = @(z) (z+1i.*u(i)).^alpha./z.^(alpha+1);
                z_values1 = linspace(beta_p, beta_p.*exp(b.*t),1000);
                z_values2 = linspace(beta_n, beta_n.*exp(b.*t),1000);
                integral_value1 = trapz(z_values1, integrand1(z_values1));
                integral_value2 = trapz(z_values2, integrand2(z_values2));
    
                values(i) = 1i.*u(i).*(1-exp(b.*t))./b.*gamma_c + ...
                    c_p.*beta_p.^alpha.*gamma(-alpha)./b.*(integral_value1 ...
                    - b.*t + alpha./beta_p.*1i.*u(i).*(1-exp(-b.*t))) + ...
                    c_n.*beta_n.^alpha.*gamma(-alpha)./b.*(integral_value2 ...
                    - b.*t - alpha./beta_n.*1i.*u(i).*(1-exp(-b.*t)));
            end
        end
    elseif strcmp(activity, 'Finite')
        aux_values = zeros(size(values));
        for i=1:numel(u)    %BISOGNA IMPLEMENTARE ANCHE LA PARTE DI C_N .38
            integrand = @(z) (z-1i.*u(i)).^alpha./z.^(alpha+1);
            z_values = linspace(beta_p, beta_p*exp(b*t));
            integral_value = trapz(z_values, integrand(z_values));
            
            aux_values(i) = 1./(b*t).*integral_value;
        end
        lamb = c_p.*beta_p*gamma(-alpha);
        values = log((exp(lamb.*t.*aux_values)-1)./(exp(lamb.*t)-1));
    end
end

end