
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

%% post process

gwParam.depthLimit = 50;
N = 30;
[s_gw, gw_M] = gen_water_growth_states(gwParam);
s_samples = cell(gw_M, N);
k_samples = cell(gw_M, N);
drawdown = cell(gw_M, N);

for i = 1:20
    filename = strcat('samples_*_', num2str(i), '.mat');
    file = dir(filename);
    load(file.name)
    index_s1 = find(s == s_gw);
    s_samples{index_s1, t} = sample_logs;
    k_samples{index_s1, t} = sample_logk;
    drawdown{index_s1, t} = dd;
end


