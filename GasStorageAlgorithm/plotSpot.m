function plotSpot(row_indices, S, type, save)
% Plot the evolution of the Spot simulation for some samples and save the
% plots for reports.
%
% INPUT:
% row_indices:      idxs of the samples to plot, to give a perfect square length
% S:                matrix with the simulations for M paths
% type:             string with the type of simualtions, just for the title
% save:             0/1 to save the plot

n = length(row_indices);

figure;

for i = 1:n
    subplot(sqrt(n), sqrt(n), i);
    plot(S(row_indices(i), :));
    grid on;
    xlabel('Days');
    ylabel('Spot Price');
    title(['Scenario ', num2str(row_indices(i)), ' ', type]);
end

% Save the plot if asked.
if save == 1
    name = sprintf('plotPrice_%s.eps', type);
    print('-depsc2', name);
end
end