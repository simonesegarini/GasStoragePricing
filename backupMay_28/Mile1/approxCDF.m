function [xgrid, CDF_hat] = approxCDF(cdf_raw, z) 
% Find the largest set of consecutive points where CDF is monotone increasing and within [0, 1]
%
% INPUT:
% cdf_raw:              raw approximation of the cdf from the FFT
% z:                    raw xgrid of the cdf
%
% OUTPUT:
% xgrid:                xgrid of the cdf
% CDF_hat:              approximated cdf


% Initialize variables to keep track of the largest set
max_length = 0;
start_index = 1;
end_index = 1;

% Iterate through cdf_raw to find the largest set
current_length = 1;
for i = 2:length(cdf_raw)
    if cdf_raw(i) >= cdf_raw(i-1) && cdf_raw(i) <= 1 && cdf_raw(i) >=0
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

% Extract the x-grid values corresponding to the largest set of consecutive points
xgrid = z(start_index:end_index);
CDF_hat = cdf_raw(start_index:end_index);

end