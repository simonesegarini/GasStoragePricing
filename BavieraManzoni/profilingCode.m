% MILESTONE 1
clear all, close all, clc

warningID = 'MATLAB:legend:IgnoringExtraEntries';
warning('off', warningID);

profile on

% SIMULATION PARAMETERS
T = 1/12;
M = 1;
N = 1e7;

% SPOT PARAMETERS
X0 = 0;

% OU-NTS VARYING ALPHA 
%alphas = [0.8, 0.6, 0.4, 0.2, -1.0, -2.0];
alphas = [0.8, 0.6, 0.4, 0.2];
b = 0.2162; sigma = 0.201; k = 0.256; theta = 0.1;

TCumulantsOUNTSalpha = zeros(numel(alphas), 4);
ECumulantsOUNTSalpha = zeros(numel(alphas), 4);

for i=1:numel(alphas)
    Xt = fgmc(X0, [alphas(i), b, sigma, k, theta], N, M, T, 'OU-NTS');
    %[TCumulantsOUNTSalpha(i,:), ECumulantsOUNTSalpha(i,:)] = computeCumulants(Xt(:,end), [alphas(i), b, sigma, k, theta], T, 'NTS-OU');
end

profile off
profile viewer