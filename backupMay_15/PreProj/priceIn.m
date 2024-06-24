function cashflows = priceIn(S, cashflows, h, N, M, delta, alpha, T, maxInjection, maxWithdraw)
% IN pricing of the storage contract with BW induction
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

% Allocate for speed efficiency
temp_cf = cashflows;

for j = T:-1:1 % bw iterations in time

    %in j i have time T+1 but S_T, if j=1 i have t = 1 but S_0
    A = [ones(length(S(:,j+1)),1), S(:,j+1), S(:,j+1).^2, S(:,j+1).^3]; 

    % OLS for CV
    b = cashflows;
    beta = A\b;
    CV = (A*beta)';

    % check all the plausible deltaV and which one has the highest exp value
    for i = 1:N
        maxAction = max(maxInjection);
        minAction = min(maxWithdraw);
        
        % Check allowed actions
        actions = maxAction/alpha:-1:minAction/alpha;
        
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
    
        for z = 1:M
            % Actions payoff + future cashflow
            total_payoff = h(S(z,j+1), actions'.*alpha) + CV(i-actions, z);
            
            % Determine the action with the highest expected future cash flow
            [~, idx] = max(total_payoff);
            
            % Update cashflow bw
            temp_cf(z,i) = exp(-delta).*cashflows(z,i-actions(idx)) + h(S(z,j+1), actions(idx).*alpha);
        end
    end

    cashflows = temp_cf;
end
end