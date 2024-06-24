function I =computeIntegral(f, interpolationGrid, numericalParams)
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

fu= f(u); % f evalued in x 

fj=fu.*exp(-1i*numericalParams.x1*numericalParams.du*(0:(N-1))); % fj 

FFT=fft(fj); % Fast Fourier Transform of fj

Curve=numericalParams.du*exp(-1i*numericalParams.u1*x).*FFT; % moltiplication of the prefactor

I=real(interp1(x,Curve,interpolationGrid));
    
end %function computeIntegral
