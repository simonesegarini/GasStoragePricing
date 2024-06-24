function values = phiXTS(us, params)
%   ADD COMMENTS


alpha = params(1); beta_p = params(3);
beta_n = params(4); c_p = params(5); c_n = params(6); gamma_c = params(7);


values = 1i.*us.*gamma_c +...
    c_p.*gamma(-alpha).*beta_p.^alpha.*((1-1i.*us./beta_p).^alpha - 1 + 1i.*us.*alpha./beta_p) + ...
    c_n.*gamma(-alpha).*beta_n.^alpha.*((1+1i.*us./beta_n).^alpha - 1 - 1i.*us.*alpha./beta_n);

end