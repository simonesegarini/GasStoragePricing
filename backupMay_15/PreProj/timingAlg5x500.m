clear all, close all, clc

means = zeros(3,5);
stds = zeros(3,5);
CIs = zeros(2,3,5);
times = zeros(3,5);
Ms = [50, 500, 5000];
seeds = 20:24;

%% STANDARD STORAGE CONTRACT DATA
% cost of injection
a1 = 0; b1 = 0;

% cost of selling
a2 = 0; b2 = 0;

% payoff
h = @(s, deltaV) -((1+a1).*s + b1).*deltaV.*(deltaV > 0) ...
                -((1-a2).*s - b2).*deltaV.*(deltaV < 0);

% time reference
ttm = 1;
delta = 0;

% volume constraints and values
Vmin = 0;
Vmax = 250000;
V0 = 100000;
VT = 100000;

% deltaV variations constraints
Imin = @(v) max(Vmin - v, -7500);
Imax = @(v) min(Vmax - v, 2500);

% penalization functions
penFunc = @(s, v) -s.*abs(v-VT).^2;

%% 1. SPOT PRICE SIMULATION 
T = 365;
S0 = 14.88;
k = 0.05; % MR parameter
mu = @(t) 0; % deterministic function of time (can be adjusted)
%mu = @(t) 0.05*t/365; 
sigma = [0.0315, 0.0945]; % low and high volatilities, MR parameter      
alpha = 2500; % width of the volume interval
N = (Vmax-Vmin)/alpha+1; % units of discretization for the volume
dV = (Vmax:-alpha:Vmin)';
index_V0 = find(dV == V0);

vol = sigma(1);

for it=1:3
    for jt=1:5
        tic
        
        M = Ms(it); % number of simulations
        
        [X, XAV] = spotSimulation(mu, vol, k, M, T+1, seeds(jt));
        S = S0*exp(X);
        SAV = S0*exp(XAV);
        
        %% 2. ASSIGN A VALUE TO THE CONTRACT AT MATURITY
        % Matrices with the cashflows
        cashflows = penFunc(S(:,end), ones(M,1)*dV');
        cashflows_AV = penFunc(SAV(:,end), ones(M,1)*dV');
        
        % Matrices for withdraw/injection limits
        maxWithdraw = Imin(dV*ones(1,T));
        maxInjection = Imax(dV*ones(1,T));
        
        %% 3. APPLY BACKWARD INDUCTION
        cashflows = priceIn(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
        cashflows_AV = priceIn(SAV, cashflows_AV, h, N, M, delta, alpha, T, maxInjection, maxWithdraw);
    
        %% 4. PRICE IN
        % Compute the final price as the mean of accumulated cash flows at t=0 across all simulations
        [means(it,jt), stds(it, jt), CIs(:, it, jt)] = normfit(0.5*(cashflows(:,index_V0) + cashflows_AV(:,index_V0)));
        
        times(it,jt) = toc;
    end
end

%% Obtain the result of a normfit on the five different runs, 
% LOAD comparisontempi.mat
meansRuns = zeros(3,1);
stdsRuns = zeros(3,1);
CIsRuns = zeros(3,2);

for i=1:3
    [meansRuns(i,1), stdsRuns(i,1), CIsRuns(i,:)] = normfit(means(i,:));
end
