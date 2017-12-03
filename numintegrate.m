
%% 
figure; 
subplot(1,2,1) 
hist(sample_logk)
subplot(1,2,2) 
hist(sample_logs)
figure;
subplot(1,3,1)
hist(drawdown_current)
subplot(1,3,2)
hist(drawdown_next)
subplot(1,3,3)
hist(dd)

%%

S_lower = 6.09E-6; 
S_upper = 2.2E-5;
K_lower = 0.0001;
K_upper = 25;

% NN info
nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
gwParam.startingHead = 337.143;


if true
    
    % Calculate norm_c
    pdf_func = str2func('unnorm_param_pdf');
    norm_c = integral2(pdf_func,log(K_lower),log(K_upper),log(S_lower),log(S_upper))


    % Use norm_c to normalized, then plot distribution
    logk = log(K_lower):.01:log(K_upper);
    logs = log(S_lower):.01:log(S_upper);
    logk_rep = repmat(logk', 1, length(logs));
    logs_rep = repmat(logs,length(logk), 1);
    norm_p = unnorm_param_pdf(logk_rep, logs_rep) / norm_c;
    
    
%     norm_p = zeros(length(logk), length(logs));
%     for i=1:length(logk)
%         for j =1:length(logs)
%             logk_rep(i,j) = logk(i);
%             logs_rep(i,j) = logs(j);
%             norm_p(i,j) =  unnorm_param_pdf(logk(i), logs(j)) / norm_c;
%         end
%     end

    % plot joint density as surface
    figure
    surf(exp(logk), exp(logs), norm_p')
    
    % save normalized density
    save('pdf', 'norm_p');

    % integrate (sum) to get marginals
    marg_s = sum(norm_p,1);
    marg_k = sum(norm_p,2);

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



% This function calcuates the unnormalized pdf for the posterior f(K,S|h(t))
% I integrate it over the full parameter space to get the normalizing
% constant





