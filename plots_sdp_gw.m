function [ ] = plots_sdp_gw(  V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, ...
    lowestCostActionIndex, sim, plotParam, s_gw, s_expand, exp_vectors, runParam, gwParam)
% Make plots from SDP and simulation results

[~,~,N] = size(T_gw_all);
[gw_M, exp_M, ~] = size(V);

%% Define population growth and demand
if plotParam.plotInitialWaterBalance
    
    population = zeros(1,N);
    population(1) = popParam.pop_initial;
    for t = 2:N
        growthScenario = popParam.growthScenario;
        growthRate = popParam.growth.(growthScenario);
        population(t) = population(t-1) * (1 + growthRate);
    end

    % For now, assume some percentage of demand per capita comes from single well
    fraction = water.demandFraction;

    %Plot initial supply - demand balance
    
    population_low = zeros(1,N);
    population_medium = zeros(1,N);
    population_high = zeros(1,N);
    population_none = zeros(1,N);
    population_low(1) = popParam.pop_initial;
    population_medium(1) = popParam.pop_initial;
    population_high(1) = popParam.pop_initial;
    population_none(1) = popParam.pop_initial;
    for t = 2:N
        growthRate = popParam.growth.low;
        population_low(t) = population_low(t-1) * (1 + growthRate);
        growthRate = popParam.growth.medium;
        population_medium(t) = population_medium(t-1) * (1 + growthRate);
        growthRate = popParam.growth.high;
        population_high(t) = population_high(t-1) * (1 + growthRate);
        growthRate = popParam.growth.none;
        population_none(t) = population_none(t-1) * (1 + growthRate);
    end
    
    gw_Minjur = ones(1,N) * gwParam.pumpingRate;
    gw_other = ones(1,N) * gwParam.otherPumpingRate;
    desal = ones(1,N) * water.desal_capacity_initial;
    desal_exp_small = ones(1,N) * water.desal_capacity_expansion.small; 
    desal_exp_large = ones(1,N) * water.desal_capacity_expansion.large; 
    waterDemand_low = demand(water, population_low, 1:N, gwParam);
    waterDemand_medium = demand(water, population_medium, 1:N, gwParam);
    waterDemand_high = demand(water, population_high, 1:N, gwParam);
    waterDemand_none = demand(water, population_none, 1:N, gwParam);
    

    f = figure;
    ax = subplot(1,2,1);
    ax.FontSize = 6;
    a1 = area(1:N, [gw_Minjur; gw_other; desal; desal_exp_large; desal_exp_small]' ./ 1E6);
    hold on;
    plot(1:N, waterDemand_low/1E6)
    plot(1:N, waterDemand_medium/1E6)
    plot(1:N, waterDemand_high/1E6)
    plot(1:N, waterDemand_none/1E6)
    legend('Minjur GW', 'Other GW', 'Desal Current', 'Desal Exp Large', 'Desal Exp Small',  'Demand Low', 'Demand Medium', 'Demand High', 'Demand None')
    legend('Location','northwest')
    legend('boxoff')
    ylabel('MCM/y')
    xlabel('Year')
    title('Water Balance: With Minjur')
    cmap = bone(5);
    for i = 1:5
        a1(i).FaceColor = cmap(i,:);
    end
    ax = subplot(1,2,2);
    ax.FontSize = 6;
    a2 = area(1:N, [gw_other; desal; desal_exp_large; desal_exp_small]' ./ 1E6);
    hold on;
    plot(1:N, waterDemand_low/1E6)
    plot(1:N, waterDemand_medium/1E6)
    plot(1:N, waterDemand_high/1E6)
    plot(1:N, waterDemand_none/1E6)
    legend('Other GW', 'Desal Current',  'Desal Exp Large', 'Desal Exp Small', 'Demand Low', 'Demand Medium', 'Demand High', 'Demand None')
    legend('Location','northwest')
    legend('boxoff')
    ylabel('MCM/y')
    xlabel('Year')
    title('Water Balance: Without Minjur')
    for i = 1:4
        a2(i).FaceColor = cmap(i+1,:);
    end
    
end



%% Visualize results: plot optimal policies

if plotParam.policyPlotsOn
    gw_step = s_gw(3) - s_gw(2);
    exp_step = s_expand(2) - s_expand(1);
    blues = colormap(cbrewer('seq', 'Blues', 6));
    oranges = colormap(cbrewer('seq', 'Oranges', 6));
    color = {blues(2,:), oranges(2,:), blues(4,:), oranges(4,:), blues(6,:), oranges(6,:), [0 0 0]};
%     fig = figure;
    times = [22 23 24 25 26 27 28 29 30];
     times = [10 11 12 13 14 15 16];
%     times = [1 2 3 4 5 6 7];
    times = [2 7 12 17 22 27];
    for t = 1:length(times)
        subplot(length(times),1,t)
        if t == 1
            patch(1,1,color{1}) % Just to make legend work, will be covered up later
            patch(1,1,color{2})
            patch(1,1,color{3})
            patch(1,1,color{4})
            patch(1,1,color{5})
            patch(1,1,color{6})
            leg = legend('No pump, no expand', 'Pump, no expand', 'No pump, small expand', 'Pump, small expand', ...
                'No pump, large expand', 'Pump, large expand');
%             legend('on');
        end
        for i = 1:gw_M
            for j = 1:exp_M
                x = [s_gw(i)-(gw_step/2) s_gw(i)-(gw_step/2) s_gw(i)+(gw_step/2) s_gw(i)+(gw_step/2)];
                if runParam.capacityDelay
                    s = 1:exp_M;
                    y = [s(j)-(exp_step/2) s(j)+(exp_step/2) s(j)+(exp_step/2) s(j)-(exp_step/2)];
                else
                    y = [s_expand(j)-(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)-(exp_step/2)];
                end
                if X1(i,j,times(t)) == 0 && X2(i,j,times(t)) == 0
                    colorThisState = color{1};
                elseif X1(i,j,times(t)) == 1 && X2(i,j,times(t)) == 0
                    colorThisState = color{2};
                elseif X1(i,j,times(t)) == 0 && X2(i,j,times(t)) == 1
                    colorThisState = color{3};
                elseif X1(i,j,times(t)) == 1 && X2(i,j,times(t)) == 1
                     colorThisState = color{4};
                elseif X1(i,j,times(t)) == 0 && X2(i,j,times(t)) == 2
                    colorThisState = color{5};
                elseif X1(i,j,times(t)) == 1 && X2(i,j,times(t)) == 2
                    colorThisState = color{6};
                elseif isnan(X1(i,j,times(t))) || isnan(X2(i,j,times(t)))
                    colorThisState = color{7};
                end
                patch(x,y,colorThisState)
                hold on  
            end
        end
    
        ax = gca;
        ax.FontSize = 6;
        ax.XTick = 0:5:s_gw(end);
        ax.XTickLabel = [];
        ax.YTickLabel = [];
        xlim([s_gw(1)-gw_step/2 s_gw(end)+gw_step/2])
        if runParam.capacityDelay
            ylim([0.5 10.5])
            ax.YTick = 1:exp_M;
        else
            ylim([s_expand(1)-exp_step/2 s_expand(end)+exp_step/2])
            ax.YTick = s_expand;
        end
        if t == length(times)
            ylabel('Added capacity state [MCM/y]')
            xlabel('Drawdown')
            if runParam.capacityDelay
                ax.YTickLabel = {'0, 0, 0', 'v, 0, 0', '2v, 0, 0', '3v, 0, 0', '0, v, 0', 'v, v, 0', '2v, v, 0', ...
                    '0, 2v, 0', 'v, 2v, 0', '0, 3v, 0', '0, 0, 3v'};
            else
                ax.YTickLabel = string(round(s_expand/1E6));
            end
            ax.XTickLabel = 0:5:s_gw(end);
        end
        title(strcat('Time step: ', num2str(times(t))))
        ax.XTickLabelRotation = 90;
    end
end



%% Plot heat maps
if plotParam.plotHeatMaps 
    hm1 = HeatMap(flipud(double(~stateInfeasible)), 'Title', 'Infeasible States (in black) from nn samples');
    addXLabel(hm1, 'time')
    addYLabel(hm1, 'drawdown')
    hm2 = HeatMap(flipud(numRelevantSamples), 'Title', 'Number of Relevant Samples: Red is high');
    addXLabel(hm2, 'time')
    addYLabel(hm2, 'drawdown')
    temp = permute(isnan(X1(:,1,:)),[1 3 2]);
    hm3 = HeatMap(flipud(double(~temp)), 'Title', 'Infeasible states (in black) from pruned tree');
end


%% Plot simulation results

[R,~] = size(sim.costOverTime);

if plotParam.simPlotsOn

if R == 1

    % Plot state evolution w/ actions
    figure;
    yyaxis left
    plot(1:N, sim.state_gw)
    hold on
    yyaxis right
    plot(1:N, sim.action_gw)
    xlabel('time')
    legend('Drawdown', 'pumping on?')

    % figure;
    % yyaxis left
    % plot(1:N, state_expand')
    % hold on
    % yyaxis right
    % plot(1:N, action_expand')
    % xlabel('time')
    % legend('Expansion state', 'Expansion decision')


    % Plot system performance
    figure
    subplot(1,3,1)
    plot(1:N, sim.failureProbOverTime)


    subplot(1,3,2)
    plot(1:N,sim.costOverTime/1E6);
    h = gca;
    % ylim([0 700])
    hold on
    bar(1:N, [sim.shortageCostOverTime./1E6; sim.expansionCostOverTime./1E6; sim.pumpingCostOverTime./1E6; sim.margDesalCostOverTime./1E6]', 'stacked');
    legend('Total cost', 'Shortage cost', 'Expansion Cost', 'Pumping Cost', 'Desal costs')
    title(strcat('Total cost [M$]: ', num2str(sum(sim.costOverTime)/1E6, '%.3E')))
    ylabel('M$')

    % subplot(1,2,2)
    % plot(1:N,shortageOverTime/1E6)
    % hold on
    % plot(1:N,demandOverTime/1E6)
    % plot(1:N,capacityOverTime/1E6)
    % bar(1:N, [minjurSupplyOverTime/1E6; othergwSupplyOverTime/1E6;  water.desal_capacity_initial*ones(1,N)/1E6;  expSupplyOverTime/1E6]', 'stacked');
    % legend('shortage', 'demand', 'capacity', 'minjur supply', 'other gw supply', 'existing desal supply', 'exp supply')
    % legend('Location', 'southwest')
    % ylabel('MCM/y');

    subplot(1,3,3)
    plot(1:N,sim.shortageOverTime/1E6)
    hold on
    plot(1:N,sim.demandOverTime/1E6)
    plot(1:N,sim.capacityOverTime/1E6)
    bar(1:N, [sim.minjurSupplyOverTime/1E6; sim.expSupplyOverTime/1E6]', 'stacked');
    legend('shortage', 'demand', 'capacity', 'minjur supply', 'exp supply')
    legend('Location', 'southwest')
    ylabel('MCM/y');
    
else
    % Plot hydrographs
    figure;
    indexLimit = find(sim.state_gw == -1 | sim.state_gw >= 100);
    hydrograph = 200 - sim.state_gw;
    hydrograph(indexLimit) = gwParam.depthLimit;
    plot(1:N, hydrograph)
    xlabel('Time [years]')
    ylabel('Drawdown [m]')
    title('Simulated Drawdown')
    
    % Plot expansion time distribution
    [r,c] = find(sim.expansionCostOverTime);
    % for every row, take the minimum column index and put NaN if none is found
    expTime = accumarray(r,c,[size(sim.expansionCostOverTime,1),1],@min,32);
    figure
    h = histogram(expTime,[0:32]);
    ax = gca;
    ax.XTick = [0.5:1:31.5];
    ax.XTickLabel = num2cell(0:31,1);
    ax.XTickLabel{end} = 'Never';
    xlabel('Expansion Year')
    ylabel('Frequency')
    title(strcat('Histogram of expansion time in ', num2str(R), ' simulations'))
    
    % Plot total shortage vs total cost
    totalCost = sum(sim.costOverTime,2);
    totalShortage = sum(sim.shortageOverTime,2);
    figure
    scatter(totalShortage,totalCost)
    title(strcat('Average Total Cost: ', num2str(sim.averageTotalCost, '%.3E')))
    
    % Bagplot!! Eventually, have different bags for learning vs no learning
    % and flexible vs no flexible
    
end

end

%% Show updated predictions over time
if plotParam.plotinfoOverTime
    
    figure;
    for k = 1:R
        numSamplesOverTime = sum(sampleIndexOverTime(:,:,k),1);
        hydrographs = cell(1,30);
        plotTimes = [1, 5, 10, 15, 20, 25];
        count = 1;
        for time = 1:30
            numSamples = numSamplesOverTime(time);
            indexSamples = find(sampleIndexOverTime(:,time,k));
            hydrographs{time} = zeros(numSamples,30);
            for i = 1:numSamples
                x = [repmat(K_samples(indexSamples(i)),[1,N]); repmat(S_samples(indexSamples(i)),[1,N]); [1:365:365*(N)]];
                tempHead = netscript(x, runParam.adjustOutput);
                hydrographs{time}(i,:) = tempHead(gwParam.wellIndex,:);
            end
            if ismember(time, plotTimes) && false
            subplot(2,3,count)
            y = [ones(numSamples,1)*200 hydrographs{time}];
            plot(0:30,y)
            count = count + 1;
            ylim([0 200])
            title(strcat('Time ', num2str(time)))
            end  
        end
        hold on
        if R >1
            subplot(R/2,2,k)
        end
        grp = [];
        y = [];
        for t=1:30
            y = [y hydrographs{t}(:,end)'];
            grp = [grp ones(1,numSamplesOverTime(t))*t]; 
        end
        boxplot(y,grp)
        hold on
        plot(xlim,[gwParam.depthLimit gwParam.depthLimit], 'k')
        ylim([0 200])
    end

end



end

