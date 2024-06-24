function volumePlot(rows, T, V0, policies, alpha, sigma, index_V0)
% Plot the evolution of the Spot simulation for some samples
%
% INPUT:
% rows:         idxs of the samples to plot, to give a perfect square length
% V0:           initial value of the volume
% policies:     3d matrix with all the policies through time
% alpha:        volume discretization parameter
% sigma:        volatility, used for plot purposes
% index_V0:     ref index for the starting volume

n = length(rows);
figure;

for j=1:n
    % iterate through policies to compute the volume step by step
    volumes = zeros(T+1,1);
    volumes(1) = V0;
    idx = index_V0;
    for i=1:T
        volumes(i+1) = volumes(i) + alpha*policies(i, idx, j);
        idx = idx - policies(i, idx, j);
    end
    
    % create the subplot
    subplot(sqrt(n), sqrt(n), j)
    plot(0:T, volumes, '-k', 'LineWidth', 4)
    xlabel('Time')
    ylabel('Volume')
    ylim([0, 250000])
    title(['Scenario ', num2str(rows(j)), ', sigma = ', num2str(sigma*100), '%']);
end
end