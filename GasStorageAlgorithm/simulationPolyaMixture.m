function increments = simulationPolyaMixture(shape, rate, a, nSim)
            
% Note that in the signature of rnbinom the probability is defines as 1 - p wrt the paper
negative_binomial = nbinrnd(shape, a,  [nSim, 1]);

increments = gamrnd(negative_binomial, a/rate, [nSim, 1]);
end