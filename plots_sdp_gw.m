function [ ] = plots_sdp_gw(  V, X1, X2, T_gw_all, cumTgw, lowestCost, ...
    lowestCostActionIndex, sim, plotParam, s_gw, s_expand, exp_vectors, runParam, gwParam, costParam, water )
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
    addpath('/Users/sarahfletcher/Documents/MATLAB/cbrewer');
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
    
%% Plot first drawdown level when build (when nocapacity)
if true
X2nocap = permute(X2(:,1,1:end-1),[1,3,2]);
X2nocap(1,:) = 0;
indexNan = isnan(X2nocap);
X2nocap(indexNan) = 0;
indexZeros = ~flipud(X2nocap == 0);
indexFirstZero =sum(cumprod(double(indexZeros),1)) + 3;
indexNan = sum(cumprod(~isnan(X2nocap)));
indexReplaceNan = indexFirstZero >= indexNan;
indexFirstZero(indexReplaceNan) = 1;
figure;
scatter(1:N, s_gw(indexFirstZero));
xlabel('Year')
ylabel('Head [m]')
title('Optimal expansion policy: Drawdown Threshold for Exapsnion over Time')
hold on
%line([0 30], [gwParam.startingHead-s_gw(end)-5 , gwParam.startingHead + 5], 'Color', 'k')
%line([0 30], [200-gwParam.depthLimit 200-gwParam.depthLimit], 'Color', 'r', 'LineStyle','--')
    
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
    yLarge = gca;
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
    indexLimit = find(sim.state_gw == -1 | sim.state_gw >= gwParam.depthLimit);
    hydrograph = gwParam.startingHead - sim.state_gw;
    hydrograph(indexLimit) = gwParam.startingHead - gwParam.depthLimit;
    plot(1:N, hydrograph)
    hold on 
    line([0 30], [gwParam.startingHead-gwParam.depthLimit gwParam.startingHead-gwParam.depthLimit ], 'Color', 'k', 'LineStyle', '--')
    xlabel('Time [years]')
    ylabel('Head [m]')
    title('Simulated Drawdown')
    ylim([0 gwParam.startingHead])
    
%     % Plot hydrograph confidence interval
%     p5 = prctile(hydrograph,5);
%     p95 = prctile(hydrograph,95);
%     med = prctile(hydrograph,50);
%     maxh = max(hydrograph);
%     minh = min(hydrograph);
%     x=1:N;
%     X=[x,fliplr(x)];                %#create continuous x value array for plotting
%     Y=[p5,fliplr(p95)];              %#create y values for out and then back
% %     fill(X,Y,'b', 'FaceAlpha', .1);figure
%     figure
%     hold on
%     plot(1:N, med, 'k');
%     Y=[minh,fliplr(maxh)];
%     fill(X,Y,'b', 'FaceAlpha', .1);
    
    % Show confidence interval evolve over time
    sample = randsample(R,1);
    head = sim.state_gw(sample,:);
    
    
    
    % Plot expansion time distribution
    [~,~, largeCost,~,~,~,~,~,~,~] = supplyAndCost( 0, 2, 0, 0, costParam, water, gwParam, 1, gwParam.pumpingRate, runParam.capacityDelay, exp_vectors);
    [~,~, smallCost,~,~,~,~,~,~,~] = supplyAndCost( 0, 1, 0, 0, costParam, water, gwParam, 1, gwParam.pumpingRate, runParam.capacityDelay, exp_vectors);
    indexLarge = sim.expansionCostOverTime == largeCost;
    indexSmall = sim.expansionCostOverTime == smallCost;
    expLargeOverTime = zeros(size(sim.expansionCostOverTime));
    expSmallOverTime = zeros(size(sim.expansionCostOverTime));
    expLargeOverTime(indexLarge) = 1;
    expSmallOverTime(indexSmall) = 1;
    [rLarge,cLarge] = find(expLargeOverTime);
    [rSmall,cSmall] = find(expSmallOverTime);
    % for every row, take the minimum column index and put NaN if none is found
    expTimeLarge = accumarray(rLarge,cLarge,[size(expLargeOverTime,1),1],@min,32);
    expTimeSmall = accumarray(rSmall,cSmall,[size(expSmallOverTime,1),1],@min,32);
    countNever = sum((expTimeLarge == 32 & expTimeSmall == 32));
    figure
    yLarge = histc(expTimeLarge,[0:32]);
    ySmall = histc(expTimeSmall,[0:32]);
    bar(0:31, [yLarge(1:end-1) ySmall(1:end-1)], 'stacked')
    hold on 
    bar(31, countNever, 'k')
    ax = gca;
    ax.XTick = [0:31];
    ax.XTickLabel = num2cell(0:31,1);
    ax.XTickLabel{end} = 'Never';
    xlim([0 32])
    xlabel('Expansion Year')
    ylabel('Frequency')
    if runParam.flexOn
        legend('Large plant', 'Small plant', 'Never')
    end
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
       
    [clrmp]=cbrewer('seq', 'Reds', N);
    numSamples = 1;

    sample = randsample(R,numSamples);
    netname = strcat('myNeuralNetworkFunction_', num2str(gwParam.nnNumber));
    netscript = str2func(netname); 
    load('T_gw_inputs_Dec4_wgaps','s_samples', 'k_samples')

    figure;
    finalHead = zeros(5,N,6);
    for k = 1:numSamples
        headSim = sim.state_gw(sample(k),:);

        p5 = zeros(N);
        p95 = zeros(N);
        samplesP5toP95 = cell(N,1);
        maxTime = N;
        for t = 1:N
            indexState = find(headSim(t) == s_gw);
            if indexState == 1
                maxTime = t-1;
                break;
            end
            K_samples = k_samples{indexState, t};
            S_samples = s_samples{indexState, t};
            headSamples = zeros(length(K_samples),N);
            for i = 1:length(K_samples)
                x = [repmat(K_samples(i),[1,N]); repmat(S_samples(i)  ,[1,N]); [365:365:365*(N)]];
                tempDrawdown = netscript(x, gwParam);
                headSamples(i,:) = gwParam.startingHead - tempDrawdown;
            end
           [headSamplesSorted, index] = sort(headSamples);
           headp5 = prctile(headSamplesSorted, 5, 1);
           headp25 = prctile(headSamplesSorted, 25, 1);
           headp50 = prctile(headSamplesSorted, 50, 1);
           headp75 = prctile(headSamplesSorted, 75, 1);
           headp95 = prctile(headSamplesSorted, 95, 1);
           p5(t,:) = headp5;
           p95(t,:) = headp95;
           finalHead(:,t,k) = [headp5(end); headp25(end); headp50(end); headp75(end); headp95(end)];

        end

        subplot(1,2,1)
        for t = 1:maxTime
            x = t:N;
            X=[x,fliplr(x)];
            scatter(t,gwParam.startingHead-headSim(t),'*', 'k')
            Y=[p5(t,t:end),fliplr(p95(t,t:end))];
            hold on
            fill(X,Y,clrmp(t,:)); 
            xlabel('Year')
            ylabel('Head [m]')
        end
        subplot(1,2,2)
        for t = 1:maxTime
            index_s1 = find(headSim(t) == s_gw);
            scatter(k_samples{index_s1,t},s_samples{index_s1,t},'o', 'k','MarkerFaceColor', clrmp(t,:))
            hold on
            xlabel('K [m/d]')
            ylabel('S')
        end
        
%         line([0 N], [200 - gwParam.depthLimit, 200 - gwParam.depthLimit], 'Color', 'r', 'LineStyle', '--')   
    end
    suptitle('Hydrograph Confidence Intervals Over Time')
    
    if false
    % Make movie for one sample
    fig = figure;
    for t = 1:maxTime
        subplot(1,2,1)
        x = t:N;
        X=[x,fliplr(x)];
        scatter(t,gwParam.startingHead-headSim(t),35, 'd', 'k', 'MarkerFaceColor', 'k')
        Y=[p5(t,t:end),fliplr(p95(t,t:end))];
        hold on
        fill(X,Y,clrmp(t,:), 'LineWidth', 1.5); 
        xlabel('Year')
        ylabel('Head [m]')
        set(gca,'linewidth',1.5)
        set(gca,'FontSize',10)
        title('Head predictions')
        subplot(1,2,2)
        scatter(samplesP5toP95{t}(1,:),samplesP5toP95{t}(2,:),'o', 'k','MarkerFaceColor', clrmp(t,:))
        hold on
        xlabel('K [m/d]')
        ylabel('S')
        set(gca,'linewidth',1.5)
        set(gca,'FontSize',10)
        set(gca,'Units','normalized')
        pos = get(gca, 'Position');
        set(gca, 'Position', pos + [0 0.05 0 -0.05])
        title('Parameters')
        xlim([0 45])
%         if t == 1
%             suptitle('Updated predictions over time (90% CI)')
%         end
        set(fig,'Position', [680 558 1000 750])
        %if t == 1 
            set(gca,'Units','normalized')
            xLabelHandle = get( gca ,'XLabel' );
            pos  = get( xLabelHandle , 'position' )
            pos1 = pos - [0 0.02 0]; 
            set( xLabelHandle , 'position' , pos1 );
        %end
        frames(t) = getframe(gcf);
    end

    % Save as movie
    myVideo = VideoWriter('infoOverTime.mpv4');
    open(myVideo)
    writeVideo(myVideo, frames);
    close(myVideo)
    
    end
    
    
% 
%     figure;
%     for i=1:6
%         figure;
%         finalH = fliplr(finalHead(:,:,i));
%         b = boxplot(rand(30));
%         h = findobj(gcf,'tag','Upper Whisker');
%         set(h,{'YData'},num2cell(finalH(end-1:end,:),1)');
%         h = findobj(gcf,'tag','Lower Whisker');
%         set(h,{'YData'},num2cell(finalH(1:2,:),1)')
%         h = findobj(gcf,'tag','Upper Adjacent Value');
%         set(h,{'YData'},num2cell(finalH([end end],:),1)')
%         h = findobj(gcf,'tag','Lower Adjacent Value');
%         set(h,{'YData'},num2cell(finalH([1 1],:),1)')
%         h = findobj(gcf,'tag','Box');
%         set(h,{'YData'},num2cell(finalH([2 4 4 2 2],:),1)')
%         h = findobj(gcf,'tag','Median');
%         set(h,{'YData'},num2cell(finalH([3 3],:),1)')
%         h = findobj(gcf,'tag','Outliers');
%         set(h,{'Visible'},{'off'})
%         %ylim([0 200])
% 
%     end
    
    
    


end

