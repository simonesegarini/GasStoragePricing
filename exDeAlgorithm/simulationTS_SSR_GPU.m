function TSrv = simulationTS_SSR_GPU(alpha, beta, theta, nSim)
% Simulation of a TS random variable on the GPU.
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

% Pre-allocate to save memory on the GPU.
TSrv = gpuArray.zeros(nSim, 1);
accepted = gpuArray.false(nSim, 1);  % Boolean array for accepted samples

% Loop until all samples are accepted
while ~all(accepted)

    temp = sum(accepted);
    disp(['SSR: generated = ', num2str(temp)])
    % Find indices of unaccepted samples
    to_generate = find(~accepted);

    % Generate random variables directly on the GPU
    num_to_generate = numel(to_generate);
    U = (rand(num_to_generate, 1, 'gpuArray') - 0.5) * pi;
    E = -log(rand(num_to_generate, 1, 'gpuArray'));
    S = (-theta * gamma(-alpha)).^(1/alpha) .* sin(alpha * U + 0.5 * pi * alpha) ...
        ./ cos(U).^(1/alpha) .* (cos((1 - alpha) * U - 0.5 * pi * alpha) ./ E) ...
        .^((1 - alpha) / alpha);

    % Generate a uniformly distributed random variable
    V = rand(num_to_generate, 1, 'gpuArray');

    % Check acceptance/rejection using logical indexing
    to_accept = V <= exp(-beta * S);
    accepted(to_generate(to_accept)) = true;
    TSrv(to_generate(to_accept)) = S(to_accept);
end

% Transfer the result back to the CPU (if needed).
TSrv = gather(TSrv);
end