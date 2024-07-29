function checkAntitheticCovariance(paths1, paths2)
    % Validate that inputs are the same size
    [nSim1, nPaths1] = size(paths1);
    [nSim2, nPaths2] = size(paths2);
    
    if nSim1 ~= nSim2 || nPaths1 ~= nPaths2
        error('The dimensions of the input matrices must match');
    end
    
    % Calculate the means of each path
    meanPaths1 = mean(paths1, 2);
    meanPaths2 = mean(paths2, 2);
    
    % Calculate the covariance between the paths
    covariance = cov(meanPaths1, meanPaths2);
    
    % Display the covariance
    fprintf('Covariance between the mean paths: %.4f\n', covariance(1, 2));
    
    % Check if the covariance is negative (indicative of antithetic variables)
    if covariance(1, 2) < 0
        disp('The covariance is negative, indicating antithetic variables were used correctly.');
    else
        disp('The covariance is non-negative, suggesting antithetic variables may not have been applied correctly.');
    end
end
