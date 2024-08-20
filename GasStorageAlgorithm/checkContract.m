function checkContract(contract)
% Check if the struct contract is correctly initialized.
%
% INPUT:
% contract:             struct with the contract

% Check for the Costs part.
if ~isfield(contract, 'Costs') ...
        || ~isfield(contract.Costs, 'TransactionCost') || isempty(contract.Costs.TransactionCost) ...
        || ~isfield(contract.Costs, 'SpreadCost') || isempty(contract.Costs.SpreadCost) ...
        || ~isfield(contract.Costs, 'TransactionProfit') || isempty(contract.Costs.TransactionProfit) ...
        || ~isfield(contract.Costs, 'SpreadProfit') || isempty(contract.Costs.SpreadProfit)
    error('Invalid costs setting, check contract.Costs struct.')
end

% Check for the Volume part.
if ~isfield(contract, 'Volume') ...
        || ~isfield(contract.Volume, 'Vmin') || isempty(contract.Volume.Vmin) ...
        || ~isfield(contract.Volume, 'Vmax') || isempty(contract.Volume.Vmax) ...
        || ~isfield(contract.Volume, 'V0') || isempty(contract.Volume.V0) ...
        || ~isfield(contract.Volume, 'VT') || isempty(contract.Volume.VT) ...
        || ~isfield(contract.Volume, 'MaxWithdraw') || isempty(contract.Volume.MaxWithdraw) ...
        || ~isfield(contract.Volume, 'MaxInjection') || isempty(contract.Volume.MaxInjection) ...
        || ~isfield(contract.Volume, 'Discretization') || isempty(contract.Volume.Discretization)
    error('Invalid volume setting, check contract.Volume struct.')
end

% Check for the Rates part.
if ~isfield(contract, 'Rates') ...
        || ~isfield(contract.Rates, 'InitialRate') || isempty(contract.Rates.InitialRate) ...
        || ~isfield(contract.Rates, 'Process') || isempty(contract.Rates.Process) ...
        || ~isfield(contract.Rates, 'Parameters') || isempty(contract.Rates.Parameters)
    error('Invalid rates setting, check contract.Rates struct.')
end

% Check for SpotInitial part.
if ~isfield(contract, 'SpotInitial') || isempty(contract.SpotInitial)
    error('Invalid setting, check contract.SpotInitial field.')
end

% Check for Simulations part.
if ~isfield(contract, 'Simulations') || isempty(contract.Simulations)
    error('Invalid setting, check contract.Simulations field.')
end

% Check for Maturity part.
if ~isfield(contract, 'Maturity') || isempty(contract.Maturity)
    error('Invalid setting, check contract.Maturity field.')
end

% Check for Steps part.
if ~isfield(contract, 'Steps') || isempty(contract.Steps)
    error('Invalid setting, check contract.Steps field.')
end

% Check for GasProcess part.
if ~isfield(contract, 'GasProcess') || isempty(contract.GasProcess)
    error('Invalid setting, check contract.GasProcess field.')
end

% Check for GasParameters part.
if ~isfield(contract, 'GasParameters') || isempty(contract.GasParameters)
    error('Invalid setting, check contract.GasParameters field.')
end

% Check for Regression part.
if ~isfield(contract, 'Regression') || isempty(contract.Regression)
    error('Invalid setting, check contract.Regression field.')
end

% Check for Order part.
if ~isfield(contract, 'Order') || isempty(contract.Order)
    error('Invalid setting, check contract.Order field.')
end

% Check for AV part.
if ~isfield(contract, 'AV') || isempty(contract.AV)
    error('Invalid setting, check contract.AV field.')
end

% Check for the Stretch part.
if ~isfield(contract, 'Stretch') || isempty(contract.Stretch)
    error('Invalid setting, check contract.Stretch field.');
end

% Check for the GPU part.
if ~isfield(contract, 'GPU') || isempty(contract.GPU)
    error('Invalid setting, check contract.GPU field.')
else
    if contract.GPU == 1
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