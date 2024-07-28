function TSrv = simulationTS_DR3(alpha, lambda, theta, nSim)
% Simulation of a TS random variable. 
% The algorithm used is the DR by Devroye 2009.
%
% INPUT:
% alpha:                first parameter of TS
% beta:                 second parameter of TS
% theta:                third parameter of TS
% nSim:                 number of variables to simulate
%
% OUTPUT:
% TSrv:                 TS simulated variables

% Consider Remark 3.

% Set up the parameters from Devroye 2009.
b = (1-alpha)/alpha; lambda_alpha = lambda.^alpha; 
gamm = lambda_alpha*alpha*(1-alpha); s_gamma = sqrt(gamm); 
s_pi = sqrt(pi); c1 = sqrt(pi/2); c2 = 2+c1; c3 = s_gamma.*c2;
xsi = (1+sqrt(2).*c3)./pi; psi = c3.*exp(-gamm.*pi.^2./8)./s_pi;

w1 = c1.*xsi./s_gamma; w2 = 2.*psi.*s_pi; w3 = xsi.*pi;

% Define the function A and B, set also the value of B(0).
S = @(x) sin(x)./x;
A = @(x) ((sin(alpha.*x).^alpha .* sin((1-alpha).*x).^(1-alpha)) ...
    ./ (sin(x))) .^ (1/(1-alpha));
B = @(x) 1./(A(x).^(1-alpha));
B0 = alpha^(-alpha) * (1-alpha)^(-(1-alpha));
BoverB0 = @(x) S(x)./((S(alpha.*x).^alpha) .* (S((1-alpha).*x).^(1-alpha)));

% Now working on the inner repeat of Devroye

X = zeros(nSim, 1);

for i=1:nSim
    % disp(['Simulazione numero: ', num2str(i)])
    accepted_outer = false;
    while ~accepted_outer
        
        accepted_inner = false;
        % Generate U with density proportional to g**.
        
        while ~accepted_inner
            
            % Generate V, W_first uniformly in [0,1].
            V = rand(1);
            W_first = rand(1);
        
            if gamm >= 1
                if V < w1/(w1+w2)
                    U = abs(randn(1))./s_gamma;
                else
                    U = pi.*(1-W_first.^2);
                end
            else % gamma < 1
                if V < w3/(w3+w2)
                    U = pi.*W_first;
                else
                    U = pi.*(1-W_first.^2);
                end
            end
            
            % Generate W uniformly in [0,1].
            zeta = sqrt(BoverB0(U));
            
            z = 1./(1 - (1+alpha.*zeta./s_gamma).^(-1./alpha));
        
            rho = pi.*exp(-lambda_alpha.*(1-1./zeta.^2))./(z+(1+c1).*s_gamma./zeta);
            d = 0;

            if ((U>=0) && (gamm >=1))
                d = d+xsi.*exp(-gamm.*U.^2./2);
            end
            if ((U>0) && (U<pi))
                d = d + psi./sqrt(pi-U);
            end
            if ((U>=0) && (U<=pi) && (gamm <1))
                d = d + xsi;
            end

            rho = rho.*d;

            W = rand(1);
            Z = rho*W;

            % Check acceptance/rejection.
            if ((U < pi) && (Z <= 1))
                accepted_inner = true;
            end
        end
    
        % Generate X with density proportional to g(x,U).
        % Set up constants.
        a = A(U); m = (b/a)^alpha * lambda_alpha;
        delta = sqrt(m*alpha/a);
        a1 = delta*c1;
        a3 = z/a;
        s = a1 + delta + a3;
    
        % Generate V_first uniformly in [0,1].
        V_first = rand(1);
        N_first = 0;
        E_first = 0;
    
        % Filter V_first to simulate differently X

        if V_first < a1/s
            N_first = randn(1);
            X(i) = m-delta.*abs(N_first);
        else
            if V_first < (a1+delta)/s
                U_first = rand(1);
                X(i) = m + delta.*U_first;
            else
                E_first = -log(rand(1));
                X(i) = m + delta + E_first.*a3;
            end
        end

        % Generate E as an exponential rv.
        E = -log(Z);
    
        % Check acceptance/rejection.
        c = a*(X(i) - m) + exp((1/alpha)*log(lambda_alpha) -b*log(m)) ...
            * ((m/X(i))^b-1);
        if X(i) < m
            c = c - N_first^2/2;
        elseif X(i) > m + delta
            c = c - E_first;
        end

        if ((X(i) >= 0) && (c<=E))
            accepted_outer = true;
        end
    end
end

TSrv = 1./(X.^b);
end