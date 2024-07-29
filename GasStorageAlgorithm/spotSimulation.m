function X = spotSimulation(model, params, N, M, T, Mfft, seed, scaling_fact, alg)
% Spot simulation based on model selected: OU returns 2*M simulations by
% using AV algorithm, OU-NTS and OU-TS return M simulations using fgmc
% algorithm
%
% INPUT:
% model:                char for model selection
% params:               vector with the parameters of the model
% N:                    number of simulations
% M:                    number of timesteps
% T:                    time horizon
% Mfft:                 parameter for FFT
% seed:                 set the seed of the simulation
% scaling_fact:         parameter to stretch the CDF and handle low alphas
% alg:                  algorithm selection (optional)
%
% OUTPUT:
% X:                    logprice from t = 0 to t = T+1

switch model
    case 'OU'
        X = spotSimulationOU(params, N, M, T, seed);

    case {'OU-NTS', 'NTS-OU'}
        X = fgmc(model, params, N, M, T, Mfft, seed, scaling_fact);

    case {'OU-TS', 'TS-OU'}
        if nargin == 9
            switch alg
                case 'fgmc'
                    X = fgmc(model, params, N, M, T, Mfft, seed, scaling_fact);
                case 'exde'
                    X = exDe(model, params, N, M, T, seed);
            end
        else
            X = fgmc(model, params, N, M, T, Mfft, seed, scaling_fact);
        end
end
end