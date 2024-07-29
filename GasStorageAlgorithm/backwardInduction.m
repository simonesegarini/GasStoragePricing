function [cashflows, policies, regressors] = backwardInduction(S, cashflows, h, N, M, rates, alpha, T, maxInjection, maxWithdraw, numKnots, method)
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
% numKnots:             number of knots for B-spline or polynomial degree + 1 for polynomial regression
% method:               'polynomial' or 'bspline' for regression type
%
% OUTPUT:
% cashflows:            matrix with the cashflows at t0
% policies:             3d matrix with all the policies through time
% regressors:           regression parameters needed for OUT pricing

% Parameters to give back, not given back in priceIn but useful for
% priceOut and for plotSpot.
policy = zeros(N, M);
policies = zeros(T, N, M);
regressors = zeros(numKnots, N, T-1);

% Run the procedure with two matrices to avoid using one 3d matrix,
% computationally more efficient and memory saving.
temp_cf = cashflows;

for j = T:-1:1 % bw iterations in time

    % In j I have time T+1 but S_T, if j=1 I have t = 1 but S_0 (idxs reference).

    % Perform regression to compute CV based on the selected method.
    switch method
        case 'polynomial'
            % Create the design matrix for polynomial regression.
            A = ones(length(S(:, j+1)), numKnots);
            for d = 1:(numKnots - 1)
                A(:,d+1) = S(:, j+1).^d;
            end

            % OLS for CV.
            b = cashflows;
            regressors(:, :, j) = A\b;
            CV = (A*regressors(:, :, j))';

        case 'bspline' 
            error('Procedure at the moment not available, trying to fix some errors.')
            
        otherwise
            error('Unknown regression method. Use "polynomial" or "bspline".')
    end

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
            temp_cf(z,i) = exp(-rates(z, j+1)).*cashflows(z,i-policy(i,z)) + h(S(z,j+1), policy(i,z)'.*alpha);
        end
        
    end
    cashflows = temp_cf;
    policies(j,:,:) = policy;
end
end