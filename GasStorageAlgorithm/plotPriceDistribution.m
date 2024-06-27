function plotPriceDistribution(prices, numberSimulations, model, extraParam, save)
% Plot the distribution of the prices of the storage contract as an histogram 
% and save for reports.
%
% INPUT:
% prices:           vector with the prices to plot
% type:             string with the type of simualtions, just for the title
% save:             0/1 to save the plot

figure;
histogram(prices, 'NumBins', numberSimulations/10)

% Switch to handle the title info based on the process.
switch model
    case {'OU-L', 'OU-H'}
        title(['Value distribution (' num2str(numberSimulations) ' paths), ', model, ', sigma = ', num2str(extraParam), '%, IN-SAMPLE']);
    otherwise
        title(['Value distribution (' num2str(numberSimulations) ' paths), ' model, ', alpha = ', num2str(extraParam), ', IN-SAMPLE'])
end
ylabel('Frequency')
xlabel('Value [Euro]')

% Save the plot if asked.
if save == 1
    name = sprintf('plotStoragePrices_%s.eps', model);
    print('-depsc2', name);
end
end