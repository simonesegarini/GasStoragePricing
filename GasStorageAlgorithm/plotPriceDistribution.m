function plotPriceDistribution(prices, numberSimulations, model, extraParam, save)
% Plot the distribution of the prices of the storage contract as a histogram 
% and save for reports.
%
% INPUT:
% prices:               vector with the prices to plot
% numberSimulations:    number of simulations of MC
% model:                flag with the type of simulations, just for the title
% extraParam:           value of alpha for the title in levy-ou and ou-levy plots
% save:                 0/1 to save the plot

figure;
histogram(prices, 'NumBins', max(10, numberSimulations / 10), 'Normalization', 'probability', 'FaceColor', [0.2 0.2 0.8]);
grid on;

% Switch to handle the title info based on the process.
switch model
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
    otherwise
        type = 'OU';
end

% Set the title based on the model
if model == 1
    title(['Value distribution (', num2str(numberSimulations), ' paths), OU'], 'FontSize', 12);
else
    title(['Value distribution (', num2str(numberSimulations), ' paths), ', type, ', \alpha = ', num2str(extraParam)], 'FontSize', 12);
end

ylabel('Frequency', 'FontSize', 10);
xlabel('Value [Euro]', 'FontSize', 10);

% Save the plot if requested
if save == 1
    name = sprintf('plotStoragePrices_%s.eps', type);
    print('-depsc2', name);
end
end
