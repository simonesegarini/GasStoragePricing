function plotSpot(row_indices, S, type, sigma)
% Plot the evolution of the Spot simulation for some samples
%
% INPUT:
% row_indices:  idxs of the samples to plot, to give a perfect square length
% S:            matrix with the simulations for M paths
% type:         string to check for MC or MCAV
% sigma:        volatility, used for plot purposes

n = length(row_indices);

figure;

for i = 1:n
    subplot(sqrt(n), sqrt(n), i);
    plot(S(row_indices(i), :));
    grid on;
    xlabel('Days');
    ylabel('Spot Price');
    if strcmp(type, 'MC')
        title(['Scenario ', num2str(row_indices(i)), ', sigma = ', num2str(sigma*100), '%']);
    elseif strcmp(type, 'MCAV')
        title(['Scenario ', num2str(row_indices(i)), ' AV, sigma = ', num2str(sigma*100), '%']);
    end
end

end