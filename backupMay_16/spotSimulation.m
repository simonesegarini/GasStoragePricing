function X = spotSimulation(model, params, N, M, T, Mfft, seed)

switch model
    case 'OU'
        X = spotSimulationOU(params, N, M, T, seed);

    case 'OU-NTS'
        X = fgmc(model, params, N, M, T, Mfft, seed);

    case 'OU-TS'
        X = fgmc(model, params, N, M, T, Mfft, seed);

end


end