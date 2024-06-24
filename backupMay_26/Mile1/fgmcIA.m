function increments = fgmcIA(U, params, M, dt, model, activity)
% FGMC method for Infinite Activity processes
%
% INPUT:
% U:                    matrix with the simulated Uniform RV
% params:               vector of paramteres for the specified model
% M:                    parameter for FFT
% dt:                   time step
% model:                model selected
% activity:             model activity, needed for the FA case
%
% OUTPUT:
% increment:            Z_deltaTj, stochastic increment for OU-Levy


% TS DO THE BILATERAL CASE
% Assign parameters and compute analyticity strip
if strcmp(model, 'OU-NTS') || strcmp(model, 'NTS-OU')
    alpha = params(1); b = params(2); sigma = params(3);
    k = params(4); theta = params(5);

    A_as = sqrt(theta^2 + (2*sigma^2*(1-alpha)) / k);
    p_n = (theta-A_as)/(sigma.^2);
    p_p = (theta+A_as)/(sigma.^2);
    
elseif strcmp(model, 'OU-TS') || strcmp(model, 'TS-OU')
    alpha = params(1); b = params(2); beta_p = params(3);
    beta_n = params(4); c_p = params(5); c_n = params(6);
    gamma_c = params(7);

    p_n = -beta_p;
    p_p = beta_n;
end

N = 2^M;
% Select shift & FFT parameters based on the decay of the process
if strcmp(activity, 'Finite') % Power decay
    % Finite activity we select N and an arbitrary du
    a = abs(0.25*max(-p_n, p_p));

    du = 0.0971; u_N = 0.5*du*(N-1); u_1 = -u_N; % Best for OU-NTS
    dx=2*pi/(N*du); x_1=-dx*(N-1)/2; x_N=-x_1;
elseif strcmp(activity, 'Infinite') % Exponential decay
    if alpha == 0  % GAMMA TS AND VG NTS CASES
        a = 0.25*max(-p_n, p_p);

        du = 0.0971; u_N = 0.5*du*(N-1); u_1 = -u_N;
        dx=2*pi/(N*du); x_1=-dx*(N-1)/2; x_N=-x_1;
    else
        a = 0.5*max(-p_n, p_p);
        
        % Infinite activity we select N and get du from the relation .33
        if strcmp(model, 'OU-NTS')
            omega = 2*alpha;
            l = 1/alpha * ((1-alpha)/k)^(1-alpha) * (sigma^2/2)^alpha * (1-exp(-2*alpha*b*dt))/(2*alpha*b);
        elseif strcmp(model, 'NTS-OU')
            omega = 2*alpha;
            l = 1/alpha * ((1-alpha)/k)^(1-alpha) * (sigma^2/2)^alpha * (1-exp(-2*alpha*b*dt));
        elseif strcmp(model, 'OU-TS')
            % CHECK THE COMPUTATION OF THESE PARAMETERS, MISSING CN
            omega = alpha;
            l = -(c_n+c_p)*gamma(-alpha)*cos(alpha*pi/2)*(1-exp(-alpha*b*dt))/(alpha*b);
        elseif strcmp(model, 'TS-OU')
            % CHECK THE COMPUTATION OF THESE PARAMETERS, MISSING CN
            if alpha == 1
                omega = 1;
                l = c_p*pi/2*(1-exp(-b*dt));
            else
                omega = alpha;
                l = -(c_p+c_n)*gamma(-alpha)*cos(alpha*pi/2)*(1-exp(-alpha*b*dt));
            end
        end
        du = (2*pi*abs(a)/(l*N^omega))^(1/(omega+1));
        u_1 = -(N-1)/2*du; u_N = -u_1;
        dx=2*pi/(N*du); x_1=-dx*(N-1)/2; x_N=-x_1;
    end
end
% NON SAREBBE DA FARE COSI MA NON CAPISCO PERCHE' VENGONO RISULTATI BELLI
% du = 0.0980; u_N = 0.5*du*(N-1); u_1 = -u_N;
% dx=2*pi/(N*du); x_1=-dx*(N-1)/2; x_N=-x_1;
Ra = a < 0;

disp(['Caso: ', model, ', alpha = ', num2str(alpha), ', du = ', num2str(du)])

numericalParams.M=M;
numericalParams.u1=u_1;
numericalParams.uN=u_N;
numericalParams.du=du;
numericalParams.x1=x_1;
numericalParams.xN=x_N;
numericalParams.dx=dx;

xgrid = -5:dx:5;

% FFT to retrieve the CDF on the xgrid
f = @(u) exp(LogCharFunc(u+1i.*a, dt, params, model, activity))./(1i.*(u+1i.*a));

I_fft= computeIntegral(f, xgrid, numericalParams);
RawCDF = Ra - exp(a.*xgrid)./(2*pi).*I_fft;

% %Plot of the CDF obtained by inverting the CF
% figure;
% plot(xgrid, RawCDF, '-k');
% title(['Raw CDF with alpha = ', num2str(alpha)])

% Restrict the xgrid by taking the largest set where the raw CDF is mototone incresing
[xgrid_hat, CDF_hat] = approxCDF(RawCDF, xgrid);

% Remove duplicate CDF_hat values and corresponding xgrid values
[CDF_hat_unique, idx] = unique(CDF_hat);
xgrid_unique = xgrid_hat(idx);

% PROVARE AD AGGIUNGERE INTERPOLATORE A MANO NEI PUNTI SEGUENDO LE
% DECADENZE COME FATTO IN FE QUEST'ANNO
first_last = [xgrid_unique(1), CDF_hat_unique(1); 
              xgrid_unique(end), CDF_hat_unique(end)];
% Invert the CDF using a spline interpolation
inverse_interpolator = @(y) interp1(CDF_hat_unique, xgrid_unique, y, 'spline');
increments = inverse_interpolator(U);

idxs_below = find(U < CDF_hat_unique(1));
idxs_over = find(U > CDF_hat_unique(end));
disp(length(idxs_below) + length(idxs_over))

% figure;
% plot(xgrid_hat, CDF_hat, '--r')
% title(['Plot with alpha =  ', num2str(alpha)])
% hold on;
% plot(increments(1:100), U(1:100), 'ok');
% plot(first_last(:,1), first_last(:,2), '*g')
% legend('Interpolated Inverted CDF', 'Interpolated values', 'Limits')

end