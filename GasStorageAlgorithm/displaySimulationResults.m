function displaySimulationResults(simulation)
% Display the results of the simulation.
%
% INPUT:
% simulation: struct containing the simulation results including ECumulants,
%             TCumulants, time, and Alphas.

% Check that necessary fields are present in the simulation struct.
requiredFields = {'Alphas', 'ECumulants', 'TCumulants', 'times'};
for i = 1:length(requiredFields)
    if ~isfield(simulation, requiredFields{i})
        error('Missing field: %s in the simulation struct.', requiredFields{i});
    end
end

% Extract relevant fields.
Alphas = simulation.Alphas;
ECumulants = simulation.ECumulants;
TCumulants = simulation.TCumulants;
times = simulation.times;

% Display the results for each Alpha.
fprintf('Simulation Results:\n');
fprintf('===================\n');

for i = 1:length(Alphas)
    fprintf('Alpha: %.1f\n', Alphas(i));
    fprintf('ECumulants: [%.3f, %.3f, %.3f, %.3f]\n', ECumulants(i, :));
    fprintf('TCumulants: [%.3f, %.3f, %.3f, %.3f]\n', TCumulants(i, :));
    fprintf('Time: %.2f\n', times(i));
    fprintf('-------------------\n');
end
end
