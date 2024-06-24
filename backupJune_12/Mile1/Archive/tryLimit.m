clear all, close all, clc

%%
% parameters
alpha = 0.4; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;
% an strip
p_n = -beta_p;
p_p = beta_n;
dt = 1/12;
a = 0.5*max(-p_n, p_p);

zs = linspace(beta_n, beta_n*exp(b*dt));
us = -500:1:1500;

%% exact int

integrand_pos = @(z, u) ((z+a-1i.*u).^alpha)./(z.^(alpha+1));
integrand_neg = @(z, u) ((z-a+1i.*u).^alpha)./(z.^(alpha+1));
values_int_pos = zeros(size(zs));
values_int_neg = zeros(size(zs));
for i=1:numel(us)
    values_int_pos(i) = trapz(zs, integrand_pos(zs, us(i)));
    values_int_neg(i) = trapz(zs, integrand_neg(zs, us(i)));
end

values = 1i.*(us+1i.*a).*(1-exp(-b*dt))*gamma_c + ...
    c_p.*beta_p.*gamma(-alpha)./b.*(values_int_pos-b*dt+alpha./beta_p.*1i.*(us+1i.*a).*(1-exp(-b*dt))) + ...
    c_n.*beta_n.*gamma(-alpha)./b.*(values_int_neg-b*dt-alpha./beta_p.*1i.*(us+1i.*a).*(1-exp(-b*dt)));

%% approx value

approx_int_pos = cos(alpha.*pi/2).*(1-exp(-alpha.*b.*dt))./(alpha.*beta_p.^alpha).*abs(us).^alpha;
approx_int_neg = cos(alpha.*pi/2).*(1-exp(-alpha.*b.*dt))./(alpha.*beta_n.^alpha).*abs(us).^alpha;

values_approx = 1i.*(us+1i.*a).*(1-exp(-b*dt))*gamma_c + ...
    c_p.*beta_p.*gamma(-alpha)./b.*(approx_int_pos-b*dt+alpha./beta_p.*1i.*(us+1i.*a).*(1-exp(-b*dt))) + ...
    c_n.*beta_n.*gamma(-alpha)./b.*(approx_int_neg-b*dt-alpha./beta_p.*1i.*(us+1i.*a).*(1-exp(-b*dt)));
%%
plot(us, real(values))
hold on
plot(us, real(values_approx))
legend('Integrated values', 'Approximated values')
