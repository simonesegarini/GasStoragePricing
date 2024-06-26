function [xgrid, CDF_hat] = approxCDF(CDF_raw, xgrid_raw, toll)
% Given the points of a function, respectively (xgrid_raw, cdf_raw), find 
% the largest set of consecutive points where cdf_raw is monotone increasing
% and within
% [0, 1]. 
% The result are given back in (xgrid, CDF_hat), the function
% serves as a filter to obtain a Cumulative Distribution Function from our
% data.
%
% INPUT:
% cdf_raw:              raw approximation of the cdf from the FFT
% xgrid_raw:            raw xgrid of the cdf
% toll:                 tollerance for CDF selection
%
% OUTPUT:
% xgrid:                xgrid of the cdf
% CDF_hat:              approximated cdf

% Initialize variables to keep track of the largest set.
max_length = 0;
start_index = 1;
end_index = 1;

% Iterate through cdf_raw to find the largest set.
current_length = 1;
for i = 2:length(CDF_raw)
    if CDF_raw(i) >= CDF_raw(i-1)-toll && CDF_raw(i) <= 1 && CDF_raw(i) >= 0
        current_length = current_length + 1;
        if current_length > max_length
            max_length = current_length;
            start_index = i - max_length + 1;
            end_index = i;
        end
    else
        current_length = 1;
    end
end

% Extract the x-grid values corresponding to the largest set of consecutive
% points.
xgrid = xgrid_raw(start_index:end_index);
CDF_hat = CDF_raw(start_index:end_index);

% Check if the first value of CDF_hat is negative.
if CDF_hat(1) < 0
    % Find the first non-negative value
    first_positive_index = find(CDF_hat >= 0, 1);
    if ~isempty(first_positive_index)
        % Adjust the xgrid and CDF_hat to start from the first non-negative
        % value.
        xgrid = xgrid(first_positive_index:end);
        CDF_hat = CDF_hat(first_positive_index:end);
    else
        % If all values are negative, return empty arrays.
        xgrid = [];
        CDF_hat = [];
    end
end
end