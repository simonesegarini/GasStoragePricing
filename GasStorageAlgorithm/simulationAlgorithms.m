function [X, ECumulants, TCumulants, times] = simulationAlgorithms(simulation)
% Simulate a process given the desired algorithm and compute the cumulants
% and check the time used.
%
% INPUT:
% simulation:           struct with the simulation info
%
% OUTPUT:
% X:                    evolution of the process
% ECumulants:           empirical cumulants
% TCumulants:           theoretical cumulants
% time:                 computational time

% Check if a seed is given.
if isfield(simulation, 'Seed') && ~isempty(simulation.Seed)
    seed = simulation.Seed;
else
    seed = randi(10000); % Generate a random seed between 1 and 10000.
end

rng(seed)

% Check general struct of the contract.
checkSimulation(simulation)

% Allocate the output variables.
TCumulants = zeros(numel(simulation.Alphas), 4);
ECumulants = zeros(numel(simulation.Alphas), 4);
times = zeros(numel(simulation.Alphas), 1);

% Check if an unfeasible combination is given.
if any(simulation.Alphas >= 1) && (simulation.Process == 22 || simulation.Process == 23 || simulation.Process == 32 || simulation.Process == 33)
    error('DR & SSR are not available for values of alpha greater or equal than 1.')
end

for i=1:numel(simulation.Alphas)
    disp(['Processing alpha = ', num2str(simulation.Alphas(i))])
    fprintf('===================\n');
    switch simulation.Process
        case 21
            tic
            X = fgmc('OU-TS', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 16, simulation.Stretch);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'OU-TS');
        case 22
            tic
            X = exDe('OU-TS', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 'SSR', simulation.GPU);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'OU-TS');
        case 23
            tic
            X = exDe('OU-TS', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 'DR', simulation.GPU);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'OU-TS');
        case 31
            tic
            X = fgmc('TS-OU', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 16, simulation.Stretch);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'TS-OU');
        case 32
            tic
            X = exDe('TS-OU', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 'SSR', simulation.GPU);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'TS-OU');
        case 33
            tic
            X = exDe('TS-OU', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 'DR', simulation.GPU);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'TS-OU');
        case 4
            tic
            X = fgmc('OU-NTS', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 16, simulation.Stretch);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'OU-NTS');

        case 5
            tic
            X = fgmc('NTS-OU', [simulation.Alphas(i), simulation.Parameters], ...
                simulation.Simulations, simulation.Steps, simulation.Maturity, 16, simulation.Stretch);
            times(i) = toc;
            [TCumulants(i,:), ECumulants(i,:)] = computeCumulants(X(:,end), [simulation.Alphas(i), simulation.Parameters], simulation.Maturity, 'NTS-OU');
    end
end