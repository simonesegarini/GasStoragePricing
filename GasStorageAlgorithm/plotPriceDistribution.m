function plotPriceDistribution(prices, numberSimulations, model, extraParam, save)
% Plot the distribution of the prices of the storage contract as an histogram 
% and save for reports.
%
% INPUT:
% prices:               vector with the prices to plot
% numberSimulations:    number of simulations of MC
% model:                flag with the type of simulations, just for the title
% extraParam:           value of alpha for the title in levy-ou and ou-levy plots
% save:                 0/1 to save the plot

figure;
histogram(prices, 'NumBins', numberSimulations/10)

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
end

switch model
    case 1
        title(['Value distribution (' num2str(numberSimulations) ' paths), OU']);
    otherwise 
        title(['Value distribution (' num2str(numberSimulations) ' paths), ' type, ', alpha = ', num2str(extraParam)])
end
ylabel('Frequency')
xlabel('Value [Euro]')

% Save the plot if asked.
if save == 1
    name = sprintf('plotStoragePrices_%s.eps', type);
    print('-depsc2', name);
end
end