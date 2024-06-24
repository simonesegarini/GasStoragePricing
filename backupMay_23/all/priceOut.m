function cashflows = priceOut(S, cashflows, betas, h, N, M, delta, alpha, T, maxInjection, maxWithdraw)
% Compute the price with the OUT method by using the regression
% coefficients computed in the IN methodology
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

temp_cf = cashflows;

for j = T:-1:1 % bw iterations in time

    %in j i have time T+1 but S_T, if j=1 i have t = 1 but S_0
    A = [ones(length(S(:,j+1)),1), S(:,j+1), S(:,j+1).^2, S(:,j+1).^3]; 
    CV = (A*betas(:,:,j))';

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
end
end