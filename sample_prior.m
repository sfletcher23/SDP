
% Samples from prior

% Parameter info
K_lower = 0.9; % [ m^2/day]
K_upper = 14;
K_mean = 4.8;
K_var = 47.15;
S_lower = 6.09E-6; 
S_upper = 2.2E-5;
L_mean = 0;
L_sigma = 5;
K_mu = log((K_mean^2)/sqrt(K_var+K_mean^2));
K_sigma = sqrt(log(K_var/(K_mean^2)+1));

% Calculate prior using lognormal dist for K, uniform for S
prior_logk_samples = normrnd(K_mu, K_sigma, 1, 5000);
prior_logs_samples = unifrnd(S_lower, S_upper, 1, 5000);


% % NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;

t = 0;
s = round(gwParam.startingHead);

% Calculate head at next period for each sample
input = [sample_logk; sample_logs; repmat(1, [1 length(sample_logk)])];
drawdown_current = netscript(input, gwParam);
input = [sample_logk; sample_logs; repmat(365*(t+1), size(sample_logk))];
drawdown_next = netscript(input, gwParam);
dd = drawdown_next - drawdown_current;

figure;
subplot(1,3,1)
hist(prior_logk_samples)
subplot(1,3,2)
hist(prior_logs_samples)
sample_logk = prior_logk_samples;
sample_logs = prior_logs_samples;
subplot(1,3,3)
hist(dd)

filename = 'samples_prior';
save(filename, 'drawdown_current', 'drawdown_next', 'dd', 'sample_logk', 'sample_logs', 's', 't')
