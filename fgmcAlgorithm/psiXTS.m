function values = psiXTS(us, params)
% Used as auxiliary function to compute the CF of a TS-OU process. The CF
% is the difference of the CF of this function with different inputs.
%
% INPUT:
% us:                   values where the process has to be evaluated
% params:               parameters of the TS process
%
% OUTPUT:
% values:               values of phi_X computed in us

% Assign parameters.
alpha = params(1); beta_p = params(3);
beta_n = params(4); c_p = params(5); c_n = params(6); gamma_c = params(7);

% Evaluate the TS process.
values = 1i.*us.*gamma_c +...
    c_p.*gamma(-alpha).*beta_p.^alpha.*((1-1i.*us./beta_p).^alpha - 1 + 1i.*us.*alpha./beta_p) + ...
    c_n.*gamma(-alpha).*beta_n.^alpha.*((1+1i.*us./beta_n).^alpha - 1 - 1i.*us.*alpha./beta_n);

end