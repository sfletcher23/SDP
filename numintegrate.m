
%% 
if false
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
end
%% post process

gwParam.depthLimit = 50;
N = 30;
[s_gw, gw_M] = gen_water_growth_states(gwParam);
s_samples = cell(gw_M, N);
k_samples = cell(gw_M, N);
drawdown = cell(gw_M, N);

files = dir('sample_data');
for i = 3:length(files)
    load(strcat('sample_data/',files(i).name))
    index_s1 = find(s == s_gw);
    s_samples{index_s1, t} = sample_logs;
    k_samples{index_s1, t} = sample_logk;
    drawdown{index_s1, t} = dd;
end

save('T_gw_inputs', 's_samples', 'k_samples', 'drawdown')

%% Calculate T_gw_row
if false
s1 = 22;
t = 1;
index_s1 = find(s1 == s_gw);
dd_input = drawdown{index_s1, t};
[T_gw_row] = gw_transrow_numint(gwParam, s1, s_gw, dd_input )

end

