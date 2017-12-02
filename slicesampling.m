
% Notes: limiting distribution on S (bc uniform over narrow range) means
% less variation in K.  

%% Generate some sample paths

N = 30;
S_lower = 6.09E-6; 
S_upper = 2.2E-5;
K_lower = 0.9;
K_upper = 14;
gwParam.startingHead = 337.143;
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 

k = K_lower;
s = S_lower;
input = [repmat(log(k), [1 N]); repmat(log(s), [1 N]); 365:365:365*N];
drawdown_max = netscript(input, gwParam);

k = K_upper;
s = S_upper;
input = [repmat(log(k), [1 N]); repmat(log(s), [1 N]); 365:365:365*N];
drawdown_min = netscript(input, gwParam);
figure
plot(1:N, gwParam.startingHead - drawdown_min)
hold on
plot(1:N,gwParam.startingHead - drawdown_max)


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
% 
% k = 3;
% s = 1E-5;
% input = [repmat(log(k), [1 N]); repmat(log(s), [1 N]); 365:365:365*N];
% drawdown3 = netscript(input, gwParam);
% plot(1:N, gwParam.startingHead - drawdown3)




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

% drawdown 3:
% t=5: 62.6876
% t=10: 85.0936
% t=15: 100.4990
% t=20: 113.3558
% t=25: 122.9495
% t=30: 132.4284

%% MCMC to estimate posterior

plotOn = false;

if false

% Slice sampling
N = 10000;
pdf_func = str2func('unnormalized_pdf');
x = slicesample([1 1E-5], N, 'pdf',pdf_func,'thin',5,'burnin',1000);
save(strcat('slice_data', getenv('SLURM_JOB_ID')), x)

prior_k = lognrnd(K_mu, K_sigma,[N 1]);
prior_s = unifrnd(S_lower, S_upper,[N 1]);

end

if plotOn
figure;
% plot K
subplot(1,2,1)
[binheight,bincenter] = hist(x(:,1),[0:.5:100]);
h = bar(bincenter,binheight,'hist');
h.FaceColor = [.8 .8 1];
h.FaceAlpha = 0.5;
hold on 
[binheight,bincenter] = hist(prior_k,[0:.5:100]);
h = bar(bincenter,binheight,'hist');
h.FaceColor = [1 .8 .8];
h.FaceAlpha = 0.5;
xlim([0 100])
% plot S
subplot(1,2,2)
[binheight,bincenter] = hist(x(:,2));
h = bar(bincenter,binheight,'hist');
h.FaceColor = [.8 .8 1];
h.FaceAlpha = 0.5;
hold on 
[binheight,bincenter] = hist(prior_s);
h = bar(bincenter,binheight,'hist');
h.FaceColor = [1 .8 .8];
h.FaceAlpha = 0.5;

end


function [p] = unnormalized_pdf(params)

k = params(1);
s = params(2);

s1 = 179;

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
% if k > 15
%     prior_k = 0;
% else
    prior_k = lognpdf(k, K_mu, K_sigma);
% end
prior_s = unifpdf(s, S_lower, S_upper);

% Log transform data for nn
k = log(k);
s = log(s);

% Calculate likelihood using model
input = [k; s; 365*10];
drawdown_t_current = netscript(input, gwParam);
y = drawdown_t_current;
u = s1; 
likelihood = normpdf(y, u, L_sigma);

% Multiply prior times likelihood
p = prior_k * prior_s * likelihood; 

end
