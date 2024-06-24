function [I1, I2] =computeIntegral(interpolationGrid, numericalParams, params, a, dt)
% Compute integral using FFT method and interpolate a grid to extract a
% curve
% INPUT: 
% f:                    integrand as a function handle
% interpolationGrid:    moneyness of interest as a vector
% numericalParams:      parameters of the numerical method as a struct
%
% OUTPUT: 
% I:                    integral values

% FFT  
N=2^numericalParams.M; % number of intervals

u=linspace(numericalParams.u1,numericalParams.uN,N); % x discretization
x=linspace(numericalParams.x1,numericalParams.xN,N); % z discretization

[f1, f2] = LogCharFunc(u + 1i .* a, dt, params);

fu_1 = exp(f1) ./ (1i .* (u + 1i .* a));
fu_2 = exp(f2) ./ (1i .* (u + 1i .* a)); % f evalued in x 

fj_1=fu_1.*exp(-1i*numericalParams.x1*numericalParams.du*(0:(N-1)));
fj_2=fu_2.*exp(-1i*numericalParams.x1*numericalParams.du*(0:(N-1))); % fj 

FFT_1=fft(fj_1);
FFT_2=fft(fj_2); % Fast Fourier Transform of fj

Curve_1=numericalParams.du*exp(-1i*numericalParams.u1*x).*FFT_1;
Curve_2=numericalParams.du*exp(-1i*numericalParams.u1*x).*FFT_2; % moltiplication of the prefactor

I1 = real(interp1(x, Curve_1, interpolationGrid));
I2 = real(interp1(x, Curve_2, interpolationGrid));
    
end %function computeIntegral

function [values, values2] = LogCharFunc(us, t, params)
alpha = params(1); b = params(2); sigma = params(3);
k = params(4); theta = params(5);

integrand = @(z, u) (0.5.*sigma.^2.*u.^2.*z.^2 - 1i*theta.*u.*z + (1-alpha)./k).^alpha./z;
z_values = linspace(exp(-b*t), 1, 1000); % Adjust the number of points as needed

[Z_mesh, U_mesh] = meshgrid(z_values, us);

evaluation = integrand(Z_mesh, U_mesh);

integral_values = trapz(z_values, evaluation, 2); % Integrate along the second dimension (columns)
values = (1 - alpha) ./ (k * alpha) .* t - (((1 - alpha) / k).^(1 - alpha)) ./ (alpha * b) .* integral_values.';
values2 = zeros(size(values));

for i=1:numel(us)
    integral_value2 = trapz(Z_mesh(i,:), evaluation(i,:));
    values2(i) = (1-alpha)./(k*alpha).*t - (((1-alpha)/k).^(1-alpha))./(alpha*b).*integral_value2;
end
end