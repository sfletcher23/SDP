%% Generate some sample paths

% figure;
% N = 30;
% k = 1;
% s = 8E-6;
% input = [repmat(log(k), [1 N]); repmat(log(s), [1 N]); 365:365:365*N];
% drawdown1 = netscript(input, gwParam);
% plot(1:N,  gwParam.startingHead - drawdown1)

% hold on
% k = 14;
% s = 2E-5;
% input = [repmat(log(k), [1 N]); repmat(log(s), [1 N]); 365:365:365*N];
% drawdown2 = netscript(input, gwParam);
% plot(1:N, gwParam.startingHead - drawdown2)

% drawdown 1:
% t=5: 132.8063
% t=10: 179.7696
% t=15: 211.3471
% t=20: 239.1236
% t=25: 264.5149
% t=30: 289.9733

% drawdown 2:
% t=5: 17.8151
% t=10: 23.9065
% t=15: 27.8651
% t=20: 29.9804
% t=25: 28.1761
% t=30: 29.9804

%%

S_lower = 6.09E-6; 
S_upper = 2.2E-5;
K_lower = 0.0001;
K_upper = 15;

% 
% pdf_func = str2func('unnorm_param_pdf');
% norm_c = integral2(pdf_func,0,15,S_lower,S_upper)
% 

if false
    
    logk = log(K_lower):.01:log(K_upper);
    logs = log(S_lower):.01:log(S_upper);
    logk_rep = repmat(logk, length(logs), 1);
    logs_rep = repmat(logs', 1, length(logk));
    p = unnorm_param_pdf(logk_rep,logs_rep);
    norm_c

    % plot joint density as surface
    figure
    surf(exp(logk), exp(logs), p/norm_c)

    % integrate (sum) to get marginals
    marg_s = sum(p,2);
    marg_k = sum(p,1);

    figure
    subplot(1,2,1)
    bincenter = logk;
    binheight = marg_k;
    h = bar(bincenter,binheight,'hist');
    h.FaceColor = [.8 .8 1];

    hold on 
    subplot(1,2,2)
    bincenter = logs;
    binheight = marg_s;
    h = bar(bincenter,binheight,'hist');
    h.FaceColor = [.8 .8 1];

end

% 
% unnorm_param_pdf(1,8E-6)/norm_c
% unnorm_param_pdf(14,2E-5)/norm_c
if true
    bins = 50:0.5:100;
    p_h = zeros(1,length(bins)); 
    for i = 1:length(bins)-1
        zmin = bins(i);
        zmax = bins(i+1);
        h_func = str2func('h_pdf');
        p_h(i) =  integral3(h_func,log(K_lower),log(K_upper),log(S_lower),log(S_upper),zmin,zmax, 'AbsTol', 1e-7, 'RelTol', 1e-5);
    end
    save(strcat('next_h_dist', getenv('SLURM_JOB_ID')),'p_h', 'bins');
    figure; 
    bincenter = bins(1:end-1) + 0.5;
    binheight = p_h(1:end-1);
    h = bar(bincenter,binheight,'hist');
end


% This function calcuates the unnormalized pdf for the posterior f(K,S|h(t))
% I integrate it over the full parameter space to get the normalizing
% constant
function [p] = unnorm_param_pdf(logk, logs)

[a, b] = size(logk);

s1 = 50;
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

% This function is used to calculate the pdf of h in the next step given
% the paramter distribution in this step. I use the previously calcualted
% posterior and the coniditional probability of the next head given the
% parameter. Integrating this gives us the marginal of h(t+1).
function [p] = h_pdf(logk, logs, h)

[a, b] = size(logk);
s1 = 50;
t = 1;
norm_c = 0.0012;

% NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;

% Calculate K,S prob using above
p_param = unnorm_param_pdf(logk,logs)/norm_c;
p_param = reshape(p_param, 1, []);

% Calculate conditional drawdowon using model
k = reshape(logk, 1, []);
s = reshape(logs, 1, []);
h = reshape(h, 1, []);
input = [k; s; repmat(365*t+1, size(k))];
dd_next = netscript(input, gwParam);
input = [k; s; repmat(365*t, size(k))];
dd_now = netscript(input, gwParam);
y = dd_next-dd_now;
% conditional = zeros(size(p_param));
% indexOne = abs(h - (s1 + y)) < 1;
% conditional(indexOne) = 1/2;
conditional = lognpdf(h - (s1 + y), 0,1);

% Multiply conditional by param
p = p_param .* conditional;
p = p_param;
p = reshape(p, a, b);

end

