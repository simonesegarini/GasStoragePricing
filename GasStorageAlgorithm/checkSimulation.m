function checkSimulation(simulation)
% Check if the struct simulation is correctly initialized.
%
% INPUT:
% simulation:           struct with the simulation parameters

% Check for the Simulations part.
if ~isfield(simulation, 'Simulations') || isempty(simulation.Simulations)
    error('Invalid setting, check simulation.Simulations field.');
end

% Check for the Maturity part.
if ~isfield(simulation, 'Maturity') || isempty(simulation.Maturity)
    error('Invalid setting, check simulation.Maturity field.');
end

% Check for the Steps part.
if ~isfield(simulation, 'Steps') || isempty(simulation.Steps)
    error('Invalid setting, check simulation.Steps field.');
end

% Check for the Seed part.
% This field is not mandatory, so only check if it exists, and if it does,
% ensure it is not empty.
if isfield(simulation, 'Seed') && isempty(simulation.Seed)
    error('Invalid setting, check simulation.Seed field.');
end

% Check for the Process part.
if ~isfield(simulation, 'Process') || isempty(simulation.Process)
    error('Invalid setting, check simulation.Process field.');
end

% Check for the Parameters part.
if ~isfield(simulation, 'Parameters') || isempty(simulation.Parameters)
    error('Invalid setting, check simulation.Parameters field.');
end

% Check for the Alphas part.
if ~isfield(simulation, 'Alphas') || isempty(simulation.Alphas)
    error('Invalid setting, check simulation.Alphas field.');
end

% Check for the Stretch part.
if ~isfield(simulation, 'Stretch') || isempty(simulation.Stretch)
    error('Invalid setting, check simulation.Stretch field.');
end

% Check for the GPU part.
if ~isfield(simulation, 'GPU') || isempty(simulation.GPU)
    error('Invalid setting, check simulation.GPU field.')
else
    if simulation.GPU == 1
        % Check if Parallel Computing Toolbox is installed
        v = ver;
        if ~any(strcmp({v.Name}, 'Parallel Computing Toolbox'))
            error('Parallel Computing Toolbox is not installed.');
        end
        
        % Check if a GPU is available
        g = gpuDeviceCount;
        if g == 0
            error('No GPU device found.');
        end
    end
end

end
