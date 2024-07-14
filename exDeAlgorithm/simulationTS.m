function TSrv = simulationTS(alpha, beta, theta, nSim)
% Simulation of the part 'X1' of the exact decomposition algorithm as a TS 
% random variable. The algorithm used is the SSR by Qu, Zhao & Dassios.
%
% INPUT:
% alpha:                first parameter of TS
% beta:                 second parameter of TS
% theta:                third parameter of TS
% nSim:                 number of variables to simulate
%
% OUTPUT:
% TSrv:                 TS simulated variables

TSrv = zeros(nSim, 1);

accepted = zeros(nSim, 1);

while sum(accepted) < nSim

    % Find idxs where we didn't accept previously.
    to_generate = find(accepted == 0);
    % temp = sum(accepted);
    % disp(['To generate = ', num2str(temp)])

    % Generate a stable rv via Zolotarev's integral representation.
    Us = (rand(size(to_generate))-0.5)*pi;
    Es = exprnd(1, size(Us));
    S = (-theta.*gamma(-alpha)).^(1/alpha) .* sin(alpha.*Us + 0.5*pi*alpha) ...
        ./ cos(Us).^(1/alpha) .* (cos((1-alpha).*Us - 0.5*pi*alpha)./Es)...
        .^((1-alpha)/alpha);
    
    % Generate a uniformly distributed rv.
    U = rand(size(Us));

    % Check aceptance/rejection
    to_accept = (U <= exp(-beta * S));

    accepted(to_generate(to_accept)) = 1;
    TSrv(to_generate(to_accept)) = S(to_accept);
end

end