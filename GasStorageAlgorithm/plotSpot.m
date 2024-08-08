function plotSpot(row_indices, S, type_FLAG, save)
% Plot the evolution of the Spot simulation for some samples and save the
% plots for reports.
%
% INPUT:
% row_indices:      idxs of the samples to plot, to give a perfect square length
% S:                matrix with the simulations for M paths
% type_FLAG:        flag with the type of simulations, just for the title
% save:             0/1 to save the plot

n = length(row_indices);

figure;
switch type_FLAG
    case 1
        type = 'OU';
    case 21 
        type = 'OU-TS FGMC';
    case 22
        type = 'OU-TS SSR';
    case 23
        type = 'OU-TS DR';
    case 31 
        type = 'TS-OU FGMC';
    case 32
        type = 'TS-OU SSR';
    case 33
        type = 'TS-OU DR';
    case 4
        type = 'OU-NTS';
    case 5
        type = 'NTS-OU';
end

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