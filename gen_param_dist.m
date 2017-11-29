function [K_samples, S_samples] = gen_param_dist(sampleSize)
% Generate N samples of the aquifer parameters in time period t given the
% information scenario number 

% Starting parameter values
K_lower = 1.2; % [ m^2/day]
K_upper = 14.3;
K_mean = 4.8;
K_var = 47.15;
S_lower = 6.09E-6;
S_upper = 2.2E-5;


% Draw samples from K (lognorm) and S (triangular) distributions
K_mu = log((K_mean^2)/sqrt(K_var+K_mean^2));
K_sigma = sqrt(log(K_var/(K_mean^2)+1));
K_samples = lognrnd(K_mu,K_sigma, [1 sampleSize]); 
S_samples = unifrnd(S_lower, S_upper, [1 sampleSize]);

end