function cashflows = priceOut(S, cashflows, betas, h, N, M, rates, alpha, T, maxInjection, maxWithdraw, numKnots, method)
% Compute the price with the OUT method by using the regression
% coefficients computed in the IN methodology.
%
% INPUT:
% S:                    simulations of Spot
% cashflows:            3d matrix with the cashflows for time, paths and vol discr
% betas:                regression parameter computed during IN pricing
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
            
            CV = (A*betas(:,:,j))';

        case 'bspline'
            error('Procedure at the moment not available, trying to fix some errors.')
            % A = spcol(1, numKnots, S(:, j+1));
            % CV = (A*betas(:, :, j))';
            
        otherwise
            error('Unknown regression method. Use "polynomial" or "bspline".')
    end

    % Check all the plausible deltaV and which one has the highest exp value.
    for i =1:N
        maxAction = max(maxInjection(:));
        minAction = min(maxWithdraw(:));

        % Create vector of actions where the values represents the quantity
        % of deltaV that can be bought/sold (ex: 2 means I can buy 2*deltaV 
        % quantity of gas, -3 I can sell 3*deltaV).
        actions = maxAction/alpha:-1:minAction/alpha;
        
        % Check the possible actions, set to zero the non allowed ones.
        % Find indices where maxInjection equals alpha*actions
        idx_inj = find(maxInjection(i,j) == alpha*actions);
        if ~isempty(idx_inj)
            actions(1:idx_inj-1) = 0;
        end
        
        % Find indices where maxWithdraw equals alpha*actions
        idx_withdraw = find(maxWithdraw(i,j) == alpha*actions);
        if ~isempty(idx_withdraw)
            actions(idx_withdraw+1:end) = 0;
        end
        
        % Iterate through the simulated path to choose the best actions and
        % compute the future cash flow.
        for z = 1:M
            % Actions payoff + future cashflow.
            total_payoff = h(S(z,j+1), actions'.*alpha) + CV(i-actions, z);
            
            % Determine the action with the highest expected future cash flow.
            [~, idx] = max(total_payoff);

            % Assign the optimal action to the policy.
            policy(i, z) = actions(idx);
            
            % Update cashflow with the backward method.
            temp_cf(z,i) = exp(-rates(z, j+1)).*cashflows(z,i-policy(i,z)) + h(S(z,j+1), policy(i,z)'.*alpha);
        end
        
    end
    % Overwrite cashflows for the cycle iterations, procedure done to avoid
    % using a 3d matrix.
    cashflows = temp_cf;
end
end