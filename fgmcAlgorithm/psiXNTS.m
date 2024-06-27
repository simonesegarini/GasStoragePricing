function values = phiXNTS(us, params)
% Used as auxiliary function to compute the CF of a NTS-OU process. The CF
% is the difference of the CF of this function with different inputs.
%
% INPUT:
% us:                   values where the process has to be evaluated
% params:               parameters of the NTS process
%
% OUTPUT:
% values:               values of phi_X computed in us

% Assign parameters.
alpha = params(1); sigma = params(3);
k = params(4); theta = params(5);

% Evaluate the NTS process.
values = (1-alpha)/(k*alpha)*(1 - (1 - 1i.*k/(1-alpha)*(theta.*us + 1i.*us.^2.*sigma^2/2)).^alpha);

end