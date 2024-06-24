function std_est = compute_std(X)
    % COMPUTE_STD Calculate the standard deviation of a vector X
    %
    % Syntax: std_est = compute_std(X)
    %
    % Input:
    %   X - A vector of numerical values
    %
    % Output:
    %   std_est - The estimated standard deviation of the vector X
    
    % Ensure X is a vector
    if ~isvector(X)
        error('Input must be a vector.');
    end
    
    % Calculate the mean of the vector
    mean_X = mean(X);
    
    % Calculate the variance (unbiased estimator)
    n = length(X);
    variance_X = sum((X - mean_X).^2) / (n - 1);
    
    % Calculate the standard deviation
    std_est = sqrt(variance_X);
end
