function plotSpot(row_indices, S, type)
% Plot the evolution of the Spot simulation for some samples
%
% INPUT:
% row_indices:  idxs of the samples to plot, to give a perfect square length
% S:            matrix with the simulations for M paths
% type:         string with the type of simualtions, just for the title

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

end