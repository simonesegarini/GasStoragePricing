clear all, close all, clc

%%

T = 1/365; % Time horizon
M = 1; % Number of steps
N = 1e7; % Number of simulations
seed = 2; % Seed for the random sampling
dt = T/M;

%%
lambdaPoisson = 1; %AKA lambda
lambdaExponential = 1; %AKA beta

k = 0.5;
a = exp(-k*dt);

shape = lambdaPoisson / k;
X = zeros(N, M+1);
polyas = zeros(N, M);
X(:,1) = 10;


for j = 1:M
    polyas(:, j) = simulationPolyaMixture(shape, lambdaExponential, a, N);
    X(:, j+1) = a*X(:, j) + polyas(:,j);
end

%%

ECUM = zeros(1, 4);
ECUM(1) = mean(polyas(:,end));
ECUM(2) = var(polyas(:,end));
ECUM(3) = skewness(polyas(:,end));
ECUM(4) = kurtosis(polyas(:,end));


TCUM = zeros(1,4);
TCUM(1)= (1-a) * shape / lambdaExponential;
TCUM(2) = (1-a^2) * shape / lambdaExponential^2;
TCUM(3) = (1-a^3) / (1-a^2)^(1.5) * 2 / sqrt(shape);
TCUM(4) = (1+a^2) / (1-a^2) * 6 / shape + 3;