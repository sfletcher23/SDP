%% Test Parameter distributions

infoScenarios = {'high_narrow','high_wide', 'medium_narrow','medium_wide','low_narrow', 'low_wide' };
infoNames = {'high narrow','high wide', 'medium narrow','medium wide','low narrow', 'low wide' };
time = 0:N;


fk = figure;
for i = 1:length(infoScenarios)
    for t = 0:N
        [K_samples,~] = gen_param_dist(infoScenarios{i}, gwParam, t, N);
        subplot(3, 2, i);
        hold on
        cdfplot(K_samples)
        title(infoNames{i})
        xlim([0 7])
    end
    if i == 6
        legend(cellstr(num2str(time(:))))
    end
end
suptitle('Hydraulic Conductivity Distributions')

fs = figure;
for i = 1:length(infoScenarios)
    for t = 0:N
        [~, S_samples] = gen_param_dist(infoScenarios{i}, gwParam, t, N);
        subplot(3, 2, i);
        hold on
        cdfplot(S_samples)
        title(infoNames{i})
        xlim([0 0.3])
    end
    if i == 6
        legend(cellstr(num2str(time(:))))
    end
end
suptitle('Storativity Distributions')