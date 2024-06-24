clear all, close all, clc

%% Parameters
sigma = 0.201; theta = 0; alpha = 1.6; k = 0.256; b = 0.1; t = 1/12; N = 16;

du = 0.0971; u_N = 0.5 * du * (N - 1); u_1 = -u_N; % Best for OU-NTS
dx = 2 * pi / (N * du); x_1 = -dx * (N - 1) / 2; x_N = -x_1;

us = linspace(u_1, u_N, 2^N);
xs = linspace(x_1, x_N, 2^N);
dx = xs(2) - xs(1);
values1 = zeros(size(us));

integrand = @(z, u) (0.5 .* sigma.^2 .* u.^2 .* z.^2 - 1i * theta .* u .* z + (1 - alpha) ./ k).^alpha ./ z;
z_values = linspace(exp(-b * t), 1, 1000); % Adjust the number of points as needed

%% METODO 1
tic 
for i = 1:numel(us)
    integral_value = trapz(z_values, integrand(z_values, us(i)));
    values1(i) = (1 - alpha) ./ (k * alpha) .* t - (((1 - alpha) / k).^(1 - alpha)) ./ (alpha * b) .* integral_value;
end
time1 = toc;
figure
plot(us, real(values1))
title('Metodo 1')

%% METODO 2 (Using Broadcasting and trapz)
tic
% Vectorized integrand evaluation using broadcasting
Z = z_values(:);  % Make sure Z is a column vector
U = us(:)';       % Make sure U is a row vector

evaluation = integrand(Z, U);

% Use trapz along the z_values dimension
integral_value2 = arrayfun(@(i) trapz(z_values, evaluation(:, i)), 1:size(evaluation, 2));

values2 = (1 - alpha) ./ (k * alpha) .* t - (((1 - alpha) / k).^(1 - alpha)) ./ (alpha * b) .* integral_value2;
time2 = toc;
figure
plot(us, real(values2))
title('Metodo 2')

% Compare results
figure
plot(us, real(values1), 'r', us, real(values2), 'b--')
legend('Metodo 1', 'Metodo 2')
title('Comparison of Metodo 1 and Metodo 2')
