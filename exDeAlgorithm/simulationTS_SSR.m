function TSrv = simulationTS_SSR(alpha, beta, theta, nSim)
% Simulation of a TS random variable. 
% The algorithm used is the SSR by Qu, Zhao & Dassios 2018.
%
% INPUT:
% alpha:                stability parameter
% beta:                 tempering parameter
% theta:                scale parameter
% nSim:                 number of variables to simulate
%
% OUTPUT:
% TSrv:                 TS simulated variables

% Pre allocate to save memory.
TSrv = zeros(nSim, 1);
accepted = zeros(nSim, 1);

while sum(accepted) < nSim

    temp = sum(accepted);
    disp(['SSR: generated = ', num2str(temp)])

    % Find idxs where we didn't accept previously.
    to_generate = find(accepted == 0);

    % Generate a stable rv via Zolotarev's integral representation.
    U = (rand(size(to_generate))-0.5)*pi;
    E = -log(rand(size(U)));
    S = (-theta.*gamma(-alpha)).^(1/alpha) .* sin(alpha.*U + 0.5*pi*alpha) ...
        ./ cos(U).^(1/alpha) .* (cos((1-alpha).*U - 0.5*pi*alpha)./E)...
        .^((1-alpha)/alpha);
    
    % Generate a uniformly distributed rv.
    V = rand(size(U));

    % Check aceptance/rejection.
    to_accept = (V <= exp(-beta * S));
    accepted(to_generate(to_accept)) = 1;
    TSrv(to_generate(to_accept)) = S(to_accept);
end

end