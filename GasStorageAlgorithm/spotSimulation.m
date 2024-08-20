function [S, SAV] = spotSimulation(S0, model, params, nSim, M, T, Mfft, GPU_flag, STRETCH, seed)
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
% GPU_flag              flag for GPU usage
% seed:                 set the seed of the simulation (optional)
%
% OUTPUT:
% X:                    logprice from t = 0 to t = T+1

if nargin == 10
    rng(seed)
end

switch model
    case 1 % OU
        Xs = spotSimulationOU(params, nSim, M, T);

    case 21 % OU-TS FGMC
        Xs = fgmc('OU-TS', params, nSim, M, T, Mfft, STRETCH);

    case 22 % OU-TS EXDE SSR
        Xs = exDe('OU-TS', params, nSim, M, T, 'SSR', GPU_flag);

    case 23 % OU-TS EXDE DR
        Xs = exDe('OU-TS', params, nSim, M, T, 'DR', GPU_flag);

    case 31 % TS-OU FGMC
        Xs = fgmc('TS-OU', params, nSim, M, T, Mfft, STRETCH);

    case 32 % TS-OU EXDE SSR
        Xs = exDe('TS-OU', params, nSim, M, T, 'SSR', GPU_flag);

    case 33 % TS-OU EXDE DR
        Xs = exDe('TS-OU', params, nSim, M, T, 'DR', GPU_flag);

    case 4 % OU-NTS FGMC
        Xs = fgmc('OU-NTS', params, nSim, M, T, Mfft, STRETCH);

    case 5 % NTS-OU FGMC
        Xs = fgmc('NTS-OU', params, nSim, M, T, Mfft, STRETCH);

    otherwise
        error('Select a valid model for underlyings simulation.')
end

X = Xs(1:nSim, :); XAV = Xs(nSim+1:end, :);
S = S0*exp(X);
SAV = S0*exp(XAV);
end