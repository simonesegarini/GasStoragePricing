function kappa = calculate_cumulants(X)
    % Calculate raw moments
    mu1 = mean(X);
    mu2 = moment(X, 2);
    mu3 = moment(X, 3);
    mu4 = moment(X, 4);

    % Calculate cumulants
    kappa1 = mu1;
    kappa2 = mu2 - mu1^2;
    kappa3 = mu3 - 3*mu2*mu1 + 2*mu1^3;
    kappa4 = mu4 - 4*mu3*mu1 - 3*mu2^2 + 12*mu2*mu1^2 - 6*mu1^4;
    kappa = [kappa1, kappa2, kappa3, kappa4]*1000;
end
