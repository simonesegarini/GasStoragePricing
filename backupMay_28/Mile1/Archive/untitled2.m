clear all; close all; clc;

% parameters
alpha = 0.4; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;
% an strip
p_n = -beta_p;
p_p = beta_n;
dt = 1/12;
a = 0.5 * max(-p_n, p_p);

zs = linspace(beta_p, beta_p * exp(b * dt));
us = -10:0.01:10;

% exact int
integrand = @(z, u) (z + a - 1i * u) ./ (z .^ (alpha + 1));
values2 = zeros(size(us));

for i = 1:numel(us)
    values2(i) = trapz(zs, integrand(zs, us(i)));
end

%disp(values2);