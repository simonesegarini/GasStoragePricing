clear all, close all, clc
%%
nSim = 1e7;
dt = 1/12;
alpha = 0.4; b = 0.1; beta_p = 2.5; beta_n = 3.5; c_p = 0.5; c_n = 1; gamma_c = 0;

a = exp(-b*dt);
B = alpha/(1-alpha);
lambda = beta_p/a;
theta = c_n*(1-a^alpha)/(alpha*b);

A = @(z) ((sin(alpha.*z).^alpha .* sin((1-alpha).*z).^(1-alpha)) ...
    ./ (sin(z))) .^ (1/(1-alpha));

%%

gAlpha = @(x) B.*x.^(-B)./pi.*integral(@(u) A(u).*exp(-A(u)./(x.^B)), 0, pi, 'ArrayValued',true);

xgrid = 0:0.0001:20;

figure
histogram(gAlpha(xgrid), 'Normalization', 'pdf')

%%

SalphaLambda = exp(lambda^alpha - lambda.*xgrid).*gAlpha(xgrid);

figure
histogram(SalphaLambda, 'Normalization', 'pdf')

%%

TS = theta.*exp(-lambda.*xgrid).*gAlpha(xgrid);

figure
histogram(TS, 'Normalization', 'pdf')
