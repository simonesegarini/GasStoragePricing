function [cashflows, policies, regressors] = backwardInduction(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw)
% Simulation the evolution of the storage and the payoff with BW induction
%
% INPUT:
% S:                    simulations of Spot
% cashflows:            3d matrix with the cashflows for time, paths and vol discr
% h:                    payoff handle function
% N:                    number of discretized volumes
% M:                    number of paths simulated
% delta:                discount factor
% alpha:                width of discretized volumes
% T:                    last day of Spot trading
% maxInjection:         matrix with max injection values for all discr volumes
% maxWithdraw:          matrix with max withdraw values for all discr volumes
%
% OUTPUT:
% cashflows:            matrix with the cashflows at t0
% policies:             3d matrix with all the policies through time
% regressors:           regression parameters needed for OUT pricing

policy = zeros(N, M);
policies = zeros(T, N, M);
regressors = zeros(4, N, T-1);

temp_cf = cashflows;

for j = T:-1:1 % bw iterations in time

    %in j i have time T+1 but S_T, if j=1 i have t = 1 but S_0
    A = [ones(length(S(:,j+1)),1), S(:,j+1), S(:,j+1).^2, S(:,j+1).^3]; 

    % OLS for CV
    b = cashflows;
    regressors(:, :, j) = A\b;
    CV = (A*regressors(:, :, j))';

    % check all the plausible deltaV and which one has the highest exp value
    for i =1:N
        maxAction = max(maxInjection(:));
        minAction = min(maxWithdraw(:));
        % check allowed actions
        actions = maxAction/alpha:-1:minAction/alpha;
        
        for act_it = 1:(maxAction/alpha + 1 - minAction/alpha)-2
            if maxInjection(i,j) == alpha*actions(act_it+1)
                actions(1:act_it) = 0;
            end

            if maxWithdraw(i,j) == alpha*actions(act_it+1)
                actions(act_it+2:end) = 0;
            end

        end

        for z = 1:M
            % actions payoff + future cashflow
            total_payoff = h(S(z,j+1), actions'.*alpha) + CV(i-actions, z);
            
            % determine the action with the highest expected future cash flow
            [~, idx] = max(total_payoff);

            % assign the optimal action to the policy
            policy(i, z) = actions(idx);
            
            % update cashflow bw
            temp_cf(z,i) = exp(-delta).*cashflows(z,i-policy(i,z)) + h(S(z,j+1), policy(i,z)'.*alpha);
        end
        
    end
    cashflows = temp_cf;
    policies(j,:,:) = policy;
end
end