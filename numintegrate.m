
%% 
if false
    
    index_s1 = 29;
    t = 25;
    
    gwParam.depthLimit = 50;
    N = 30;
    [s_gw, gw_M] = gen_water_growth_states(gwParam);
    
    sample_logk = k_samples{index_s1,t};
    sample_logs = s_samples{index_s1,t};
    dd = drawdown{index_s1, t};
    figure; 
    subplot(1,3,1) 
    hist(sample_logk)
    xlabel('logk')
    title(strcat('mean: ', num2str(mean(sample_logk)), ' sd: ',  num2str(std(sample_logk))))
    subplot(1,3,2) 
    hist(sample_logs)
    xlabel('logs')
    title(strcat('mean: ',  num2str(mean(sample_logs)), ' sd: ',  num2str(std(sample_logs))))
    subplot(1,3,3)
    hist(dd)
    title(strcat('mean: ',  num2str(mean(dd)), ' sd: ',  num2str(std(dd))))
    xlabel('dd next step')
    suptitle(strcat('s1: ', num2str(s_gw(index_s1)), ' t: ', num2str(t)));
end
%% post process

gwParam.depthLimit = 50;
N = 30;
[s_gw, gw_M] = gen_water_growth_states(gwParam);
s_samples = cell(gw_M, N);
k_samples = cell(gw_M, N);
drawdown = cell(gw_M, N);
s_gw

files = dir('sample_data');
for i = 3:length(files)
    load(strcat('sample_data/',files(i).name))
    s
    t
    if ~isnumeric(s)
        s = str2num(s)
    end
    if ~isnumeric(t)
        t = str2num(t)
    end
    index_s1 = find(s == s_gw);
    s_samples{index_s1, t} = sample_logs;
    k_samples{index_s1, t} = sample_logk;
    drawdown{index_s1, t} = dd;
end

save('T_gw_inputs_Dec4_wgaps', 's_samples', 'k_samples', 'drawdown')

%% Plot some series



%% Calculate T_gw_row
if false
s1 = 22;
t = 1;
index_s1 = find(s1 == s_gw);
dd_input = drawdown{index_s1, t};
[T_gw_row] = gw_transrow_numint(gwParam, s1, s_gw, dd_input )

end

