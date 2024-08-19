function TSrv = simulationTS_DR_GPU(alpha, lambda, nSim)
% Simulation of a TS random variable using the DR algorithm by Devroye (2009)
% with GPU acceleration and without batch processing.
%
% INPUT:
% alpha:                Stability parameter
% lambda:               Tempering parameter
% nSim:                 Number of variables to simulate
%
% OUTPUT:
% TSrv:                 TS simulated variables

% Set up the parameters from Devroye 2009.
gamma = lambda^alpha * alpha * (1-alpha);
xsi = 1/pi * ((2 + sqrt(pi/2)) * sqrt(2*gamma) + 1);
psi = 1/pi * exp(-(gamma*pi^2)/8) * (2 + sqrt(pi/2)) * sqrt(gamma*pi);
w1 = xsi * sqrt(pi/(2*gamma));
w2 = 2 * psi * sqrt(pi);
w3 = xsi * pi;
b = (1-alpha)/alpha;

accepted_outer = gpuArray.false(nSim, 1);

% Define the function A and B, set also the value of B(0).
S = @(x) sin(x)./x;
A = @(x) ((sin(alpha.*x).^alpha .* sin((1-alpha).*x).^(1-alpha)) ...
    ./ (sin(x))) .^ (1/(1-alpha));
BoverB0 = @(x) S(x)./((S(alpha.*x).^alpha) .* (S((1-alpha).*x).^(1-alpha)));

% Allocate for memory efficiency.
X = zeros(nSim, 1);

while ~all(accepted_outer)

    % Allocate for memory efficiency.
    X_temp = gpuArray.zeros(nSim - sum(accepted_outer), 1);

    % Generate U with density proportional to g**.
    accepted_inner = gpuArray.false(nSim - sum(accepted_outer), 1);
    U = gpuArray.zeros(nSim - sum(accepted_outer), 1);
    z = gpuArray.zeros(nSim - sum(accepted_outer), 1);

    % These needed jst for debugging.
    W = gpuArray.zeros(nSim - sum(accepted_outer), 1);
    rho = gpuArray.zeros(nSim - sum(accepted_outer), 1);
    
    while ~all(accepted_inner)
        
        U_temp = gpuArray.zeros(nSim - sum(accepted_inner) - sum(accepted_outer),1);

        % Find idxs where we didn't accept previously.
        to_generate_inner = find(accepted_inner == 0);
        
        % Generate V, W_first uniformly in [0,1].
        V = rand(size(to_generate_inner), 'gpuArray');
        W_first = rand(size(to_generate_inner), 'gpuArray');
    
        if gamma >= 1
            V_inf_idxs = V < w1/(w1+w2);
            V_sup_idxs = V >= w1/(w1+w2);
            
            % Can improve here
            N = randn(size(V_inf_idxs), 'gpuArray');
            U_temp(V_inf_idxs) = abs(N(V_inf_idxs))./sqrt(gamma);
            U_temp(V_sup_idxs) = pi.*(1-W_first(V_sup_idxs).^2);
        else % gamma < 1
            V_inf_idxs = V < w3/(w3+w2);
            V_sup_idxs = V >= w3/(w3+w2);
    
            U_temp(V_inf_idxs) = pi.*W_first(V_inf_idxs);
            U_temp(V_sup_idxs) = pi.*(1-W_first(V_sup_idxs).^2);
        end
        
        % Generate W uniformly in [0,1].
        W_temp = rand(size(to_generate_inner), 'gpuArray');
        zeta = sqrt(complex(BoverB0(U_temp)));
        phi = (sqrt(gamma) + alpha.*zeta).^(1/alpha);
        z_temp = phi./(phi - (sqrt(gamma)).^(1/alpha));
    
        rho_temp = ( pi .* exp(-lambda.^alpha .* (1-zeta.^(-2))) .* ...
            (xsi .* exp(-0.5.*gamma.*U_temp.^2).*(U_temp >= 0).*(gamma >= 1) + ...
            psi ./ (sqrt(complex(pi-U_temp))).*(U_temp > 0).*(U_temp < pi) + ...
            xsi.*(gamma < 1).*(U_temp <= pi).*( U_temp>= 0)) )...
            ./((1 + sqrt(pi/2)) .* sqrt(gamma)./zeta + z_temp);
    
        % Check acceptance/rejection.
        to_accept_inner = (U_temp <= pi).* (W_temp.*rho_temp <= 1);
        to_accept_inner = logical(to_accept_inner);
    
        accepted_inner(to_generate_inner(to_accept_inner)) = 1;
        U(to_generate_inner(to_accept_inner)) = U_temp(to_accept_inner);
        z(to_generate_inner(to_accept_inner)) = z_temp(to_accept_inner);
        W(to_generate_inner(to_accept_inner)) = W_temp(to_accept_inner);
        rho(to_generate_inner(to_accept_inner)) = rho_temp(to_accept_inner);
    
    end

    % Define Z, useful for later.
    Z = W.*rho;

    % Find idxs where we didn't accept previously.
    to_generate_outer = find(accepted_outer == 0);

    % Generate X with density proportional to g(x,U).
    % Set up constants.
    a = A(U); m = (b.*lambda./a).^alpha;
    delta = sqrt(m.*alpha./a);
    a1 = delta.*sqrt(pi/2);
    a2 = delta;
    a3 = z./a;
    s = a1 + a2 + a3;

    % Generate V_first uniformly in [0,1].
    V_first = rand(size(to_generate_outer), 'gpuArray');

    % Filter V_first to simulate differently X.
    V_first_i1 = V_first < a1./s;
    V_first_i2 = (V_first < (a1+a2)./s).*(V_first >= a1./s);
    V_first_i2 = logical(V_first_i2);
    V_first_i3 = V_first >= (a1+a2)./s;
    
    N_first = randn(size(V_first_i1), 'gpuArray');
    E_first = -log(rand(size(V_first_i3), 'gpuArray'));

    X_temp(V_first_i1) = m(V_first_i1) - delta(V_first_i1).*abs(N_first(V_first_i1));
    U_third = rand(size(V_first_i2), 'gpuArray');
    X_temp(V_first_i2) = delta(V_first_i2).*U_third(V_first_i2) + m(V_first_i2);
    X_temp(V_first_i3) = m(V_first_i3) + delta(V_first_i3) + E_first(V_first_i3).*a3(V_first_i3);

    % Generate E as an exponential rv.
    E = -log(Z);

    % Check acceptance/rejection.
    to_accept_outer = (X_temp >= 0).* ...
        ((a.*(X_temp - m) + lambda.*(power(complex(X_temp), complex(-b, 0)) - m.^(-b))...
        -0.5.*N_first.^2.*(X_temp < m) - E_first.*(X_temp > m + delta)) <= E);
    to_accept_outer = logical(to_accept_outer);

    accepted_outer(to_generate_outer(to_accept_outer)) = 1;
    X(to_generate_outer(to_accept_outer)) = X_temp(to_accept_outer);

end

TSrv = 1./(X.^b);
end
