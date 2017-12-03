
% Define parameter space
S_lower = 6.09E-6; 
S_upper = 2.2E-5;
K_lower = 0.9;
K_upper = 14;

% Integrate posterior to get normalizing constant
pdf_func = str2func('unnorm_param_pdf');
norm_c = integral2(pdf_func,log(K_lower),log(K_upper),log(S_lower),log(S_upper));

% Use norm_c to get pdf values of normalized posterior
logk = log(K_lower):.01:log(K_upper);
logs = log(S_lower):.01:log(S_upper);
logk_rep = repmat(logk', 1, length(logs));
logs_rep = repmat(logs,length(logk), 1);
norm_p = unnorm_param_pdf(logk_rep, logs_rep) / norm_c;


% Save
taskID = getenv('SLURM_ARRAY_TASK_ID');
jobID = getenv('SLURM_JOB_ID');
filename = strcat('samples_', num2str(jobID), '_',  num2str(taskID), '.mat');
save(filename, 'norm_p');

% This function calcuates the unnormalized pdf for the posterior f(K,S|h(t))
% I integrate it over the full parameter space to get the normalizing
% constant

function [p] = unnorm_param_pdf(logk, logs)

[a, b] = size(logk);

s1 = 7.0;
t = 1;

% NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;

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
% pd = makedist('Lognormal','mu',K_mu,'sigma',K_sigma);
% pd = truncate(pd,0,40);
% prior_k = pdf(pd, k);
if exp(logk) > 15
    prior_k = 0;
else
    prior_k = normpdf(logk, K_mu, K_sigma);
end
prior_s = unifpdf(exp(logs), S_lower, S_upper);

% reshape
k = reshape(logk, 1, []);
s = reshape(logs, 1, []);
prior_s = reshape(prior_s, 1, []);
prior_k = reshape(prior_k, 1, []);

% Calculate likelihood using model
input = [k; s; repmat(365*t, size(k))];
drawdown_t_current = netscript(input, gwParam);
y = drawdown_t_current;
u = s1; 
likelihood = normpdf(y, u, L_sigma);

% Multiply prior times likelihood
p = prior_k .* prior_s .* likelihood; 
p = reshape(p, a, b);

end