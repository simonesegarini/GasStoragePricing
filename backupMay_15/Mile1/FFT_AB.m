function CDF = FFT_AB(f, interpolationGrid, numericalParams, Ra, a)

N=2^numericalParams.M; % number of intervals

u=linspace(numericalParams.u1,numericalParams.uN,N); % x discretization
x=linspace(numericalParams.x1,numericalParams.xN,N);
CDF = zeros(size(x));

% numericalParams.M=M;
% numericalParams.u1=u_1;
% numericalParams.uN=u_N;
% numericalParams.du=du;
% numericalParams.x1=x_1;
% numericalParams.xN=x_N;
% numericalParams.dx=dx;

for xit = 1:numel(x)

    fu = exp(f(0:(N-1)+0.5).*numericalParams.du + 1i*a)./(1i.*(0:(N-1)+0.5).*numericalParams.du-a)...
        .*exp(-1i.*xit.*(0:(N-1)+0.5).*numericalParams.du);

    CDF(xit) = Ra - numericalParams.du*exp(a.*interpolationGrid)/pi * real(fu);
end