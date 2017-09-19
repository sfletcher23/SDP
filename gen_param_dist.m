function [K_samples, S_samples] = gen_param_dist(infoScenario, sampleSize, t, N)
% Generate N samples of the aquifer parameters in time period t given the
% information scenario number 

% Starting parameter values
K_lower = 0.4; % [ m^2/day]
K_upper = 2.5;
K_mean = 1.170;
K_sigma = 0.56;
S_lower = 0.02;
S_upper = 3e-1;
S_mean = 0.13;


% Define mean and stddev for each distribution over time
switch infoScenario
    case 'high_narrow'
        K_start = K_mean;
        K_end = K_upper * .8 + K_mean * .2;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma / 3;
        S_start = S_mean;
        S_end = S_upper * .8 + S_mean * .2;
        S_lower_end = S_upper * .5 + S_mean * .5;
        S_upper_end = S_upper;
        
    case 'low_narrow'
        K_start = K_mean;
        K_end = K_lower * .8 + K_mean * .2;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma / 3;
        S_start = S_mean;
        S_end = S_lower * .8 + S_mean * .2;
        S_lower_end = S_lower;
        S_upper_end = S_lower * .5 + S_mean *.5;
        
    case 'medium_narrow'
        K_start = K_mean;
        K_end = K_mean;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma / 3;
        S_start = S_mean;
        S_end = S_mean;
        S_lower_end = S_lower * .5 + S_mean *.5;
        S_upper_end = S_upper *.5 + S_mean * .5;
        
    case 'high_wide'
        K_start = K_mean;
        K_end = K_upper * .8 + K_mean * .2;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma * .9;
        S_start = S_mean;
        S_end = S_upper * .8 + S_mean * .2;
        S_upper_end = S_upper;
        S_lower_end = S_lower * .5 + S_mean *.5;
        
    case 'low_wide'
        K_start = K_mean;
        K_end = K_lower * .8 + K_mean * .2;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma * .9;
        S_start = S_mean;
        S_end = S_lower * .8 + S_mean * .2;
        S_lower_end = S_lower;
        S_upper_end = S_upper * .5 + S_mean * .5;
        
    case 'medium_wide'
        K_start = K_mean;
        K_end = K_mean;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma * .9;
        S_start = S_mean;
        S_end = S_mean;
        S_lower_end = S_lower * .9 + S_mean * .1;
        S_upper_end = S_upper * .9 + S_mean * .1;
        
    case 'full_range'
        K_start = K_mean;
        K_end = K_mean;
        K_sigma_start = K_sigma;
        K_sigma_end = K_sigma;
        S_start = S_mean;
        S_end = S_mean;
        S_lower_end = S_lower;
        S_upper_end = S_upper;
        
    otherwise
        error('invalid information scenario name')
end

% Calulate K Param values for time period t
K_step = (K_end - K_start) / N;
K_t = K_start + K_step * t;
K_sigma_step = (K_sigma_end - K_sigma_start) / N;
K_sigma_t = K_sigma_start + K_sigma_step * t;

% Calculate S Param values for time period t
S_step = (S_end - S_start) / N;
S_t  = S_start + S_step * t;
S_lower_step = (S_lower_end - S_lower) / N;
S_upper_step = (S_upper_end - S_upper) / N;
S_lower_t = S_lower + S_lower_step * t;
S_upper_t = S_upper + S_upper_step * t;

% Draw samples from K (lognorm) and S (triangular) distributions
K_mu_t = log(K_t) - (K_sigma_t ^ 2) / 2;
K_samples = lognrnd(K_mu_t,K_sigma_t, [1 sampleSize]); 
pd = makedist('Triangular','a',S_lower_t,'b',S_t,'c',S_upper_t);
S_samples = random(pd, [1 sampleSize]);


end