%% updating over time plots 


% Setup
cd('/Users/sarahfletcher/Documents/MATLAB/Repository_SDP')
addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/Repository_SDP'))
addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/figure_tools'))
load('Results/Base/Dec5_base.mat')

N = 30;
[clrmp]=cbrewer('seq', 'Reds', N);
numSamples = 1;

R = 500;
sample = randsample(R,numSamples);
sample = 158;
netname = strcat('myNeuralNetworkFunction_', num2str(gwParam.nnNumber));
netscript = str2func(netname); 

% Get K and S samples
load('T_gw_inputs_Dec4','s_samples', 'k_samples'); 
step = 2; % display for every other year to reduce visual clutter

f = figure;
FigHandle = f;
figure_width = 8;
figure_height = 6.5;
font_size = 8;
line_width = 1;
export_ppi = 600;
print_png = false;
print_pdf = false;
savename = 'Results/Base/Learning_w_shortage_double';
printsetup(FigHandle, figure_width, figure_height, font_size, line_width, export_ppi, print_png, print_pdf, savename)



finalHead = zeros(5,N,6);
for k = 1:numSamples
    headSim = sim.state_gw(sample(k),:);

    p5 = zeros(N);
    p95 = zeros(N);
    samplesP5toP95 = cell(N,1);
    maxTime = N;
    volShort = cell(N,1);
    for t = 1:step:N
        % get current state
        indexState = find(headSim(t) == s_gw);
        if indexState == 1
            maxTime = t-1;
            break;
        end
        % get samples for current state
        K_samples = k_samples{indexState, t};
        S_samples = s_samples{indexState, t};
        headSamples = zeros(length(K_samples),N);
        
        % Get head for current samples
        for i = 1:length(K_samples)
            x = [repmat(K_samples(i),[1,N]); repmat(S_samples(i)  ,[1,N]); [365:365:365*(N)]];
            tempDrawdown = netscript(x, gwParam);
            headSamples(i,:) = gwParam.startingHead - tempDrawdown;
        end 
        % Sort head samples to get CIs
       [headSamplesSorted, index] = sort(headSamples);
       headp5 = prctile(headSamplesSorted, 5, 1);
       headp25 = prctile(headSamplesSorted, 25, 1);
       headp50 = prctile(headSamplesSorted, 50, 1);
       headp75 = prctile(headSamplesSorted, 75, 1);
       headp95 = prctile(headSamplesSorted, 95, 1);
       p5(t,:) = headp5;
       p95(t,:) = headp95;
       finalHead(:,t,k) = [headp5(end); headp25(end); headp50(end); headp75(end); headp95(end)];
       
       % Use head samples to calculate shortage CIs
       shortageSamples = zeros(length(K_samples),N);
       indexShort = headSamples < (gwParam.startingHead - gwParam.depthLimit);
       yearsShort = sum(indexShort,2);
       volShort{t} = yearsShort * gwParam.pumpingRate / 1E9;
    end

    subplot(2,3,2)
    for t = 1:step:maxTime
        x = t:N;
        X=[x,fliplr(x)];
        Y=[p5(t,t:end),fliplr(p95(t,t:end))];
        hold on
        fill(X,Y,clrmp(t,:)); 
        line([0 N], [gwParam.startingHead - gwParam.depthLimit, gwParam.startingHead - gwParam.depthLimit], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2) 
        scatter(t,gwParam.startingHead-headSim(t),'h', 'k', 'MarkerFaceColor', 'k', 'SizeData', 40)
        xlabel('Year')
        ylabel('Head [m]')
        title('Hydraulic Head 90% CI') 
        legend('Head CI', 'depth limit', 'simulated obs')
        legend('boxoff')
        xticks(0:4:30)
    end
    subplot(2,3,1)
    for t = 1:step:maxTime
        index_s1 = find(headSim(t) == s_gw);
        scatter(k_samples{index_s1,t},s_samples{index_s1,t},'o', 'k','MarkerFaceColor', clrmp(t,:))
        hold on
        xlabel('log K [m/d]')
        ylabel('log S')
        ylim([-12.05 -10.7])
        xlim([0 2.75])
        title('Parameter samples')
    end
    
    subplot(2,3,3)
    shortData = zeros(5000, 15);
    count = 1;
    for t = 1:step:maxTime
        shortData(:,count) = volShort{t};
        count = count +1;
    end
    ind = maxTime:-step:1;
    boxplot(shortData,'PlotStyle','traditional','BoxStyle', 'outline', 'Symbol', '', 'Colors', 'k');
    a = get(get(gca,'children'),'children');
    boxes =  a(15*2+1:15*3);
    medians = a(16:30);
%     set(medians, 'Marker', 'o')
    for j=1:length(boxes)
        hold on
        p(j) = patch(get(boxes(j),'XData'),get(boxes(j),'YData'),clrmp(ind(j),:));
        m(j) = scatter(mean(get(medians(j),'XData')), mean(get(medians(j),'YData')), 'k', 'o', 'MarkerFaceColor', 'k', 'SizeData', 16);
    end
    uistack(p,'bottom')
    whiskers = a(15*5+1:end);
    set(whiskers, 'LineStyle', '-')
    xticks(0:2:15)
    xticklabels({'0', '4', '8', '12', '16', '20', '24', '28'})
    ylabel('30-year shortages [BCM]')
    xlabel('year')
    title('Water shortages (no desal)')
    box off
end


sample = 7;
finalHead = zeros(5,N,6);
for k = 1:numSamples
    headSim = sim.state_gw(sample(k),:);
    p5 = zeros(N);
    p95 = zeros(N);
    samplesP5toP95 = cell(N,1);
    maxTime = N;
    volShort = cell(N,1);
    for t = 1:step:N
        % get current state
        indexState = find(headSim(t) == s_gw);
        if indexState == 1
            maxTime = t-1;
            break;
        end
        % get samples for current state
        K_samples = k_samples{indexState, t};
        S_samples = s_samples{indexState, t};
        headSamples = zeros(length(K_samples),N);
        
        % Get head for current samples
        for i = 1:length(K_samples)
            x = [repmat(K_samples(i),[1,N]); repmat(S_samples(i)  ,[1,N]); [365:365:365*(N)]];
            tempDrawdown = netscript(x, gwParam);
            headSamples(i,:) = gwParam.startingHead - tempDrawdown;
        end 
        % Sort head samples to get CIs
       [headSamplesSorted, index] = sort(headSamples);
       headp5 = prctile(headSamplesSorted, 5, 1);
       headp25 = prctile(headSamplesSorted, 25, 1);
       headp50 = prctile(headSamplesSorted, 50, 1);
       headp75 = prctile(headSamplesSorted, 75, 1);
       headp95 = prctile(headSamplesSorted, 95, 1);
       p5(t,:) = headp5;
       p95(t,:) = headp95;
       finalHead(:,t,k) = [headp5(end); headp25(end); headp50(end); headp75(end); headp95(end)];
       
       % Use head samples to calculate shortage CIs
       indexShort = headSamples < (gwParam.startingHead - gwParam.depthLimit);
       yearsShort = sum(indexShort,2);
       volShort{t} = yearsShort * gwParam.pumpingRate / 1E9;
             
    end
    
    for t = maxTime+1:step:N
        yearsShort = ones(size(yearsShort)) * (N - maxTime + 1);
        volShort{t} = yearsShort * gwParam.pumpingRate / 1E9;
    end

    subplot(2,3,5)
    for t = 1:step:maxTime
        x = t:N;
        X=[x,fliplr(x)];
        Y=[p5(t,t:end),fliplr(p95(t,t:end))];
        hold on
        fill(X,Y,clrmp(t,:)); 
        line([0 N], [gwParam.startingHead - gwParam.depthLimit, gwParam.startingHead - gwParam.depthLimit], 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2) 
        scatter(t,gwParam.startingHead-headSim(t),'h', 'k', 'MarkerFaceColor', 'k', 'SizeData', 40)
        xlabel('Year')
        ylabel('Head [m]')
        title('Hydraulic Head 90% CI') 
        legend('Head CI', 'depth limit', 'simulated obs')
        legend('boxoff')
        xticks(0:4:30)
    end
    subplot(2,3,4)
    for t = 1:step:maxTime
        index_s1 = find(headSim(t) == s_gw);
        scatter(k_samples{index_s1,t},s_samples{index_s1,t},'o', 'k','MarkerFaceColor', clrmp(t,:))
        hold on
        xlabel('log K [m/d]')
        ylabel('log S')
        ylim([-12.05 -10.7])
        xlim([0 2.75])
        title('Parameter samples')
    end
    
    subplot(2,3,6)
    shortData = zeros(5000, 15);
    count = 1;
    for t = 1:step:N
        shortData(:,count) = volShort{t};
        count = count +1;
    end
    ind = N:-step:1;
    boxplot(shortData,'PlotStyle','traditional','BoxStyle', 'outline', 'Symbol', '', 'Colors', 'k');
    a = get(get(gca,'children'),'children');
    boxes =  a(15*2+1:15*3);
    medians = a(16:30);
%     set(medians, 'Marker', 'o')
    for j=1:length(boxes)
        hold on
        p(j) = patch(get(boxes(j),'XData'),get(boxes(j),'YData'),clrmp(ind(j),:));
        m(j) = scatter(mean(get(medians(j),'XData')), mean(get(medians(j),'YData')), 'k', 'o', 'MarkerFaceColor', 'k', 'SizeData', 16);
    end
    uistack(p,'bottom')
    whiskers = a(15*5+1:end);
    set(whiskers, 'LineStyle', '-')
    xticks(0:2:15)
    xticklabels({'0', '4', '8', '12', '16', '20', '24', '28'})
    ylabel('30-year shortages [BCM]')
    xlabel('year')
    title('Water shortages (no desal)')
    box off
end



suptitle('Bayesian Learning Over Time')
f.PaperSize = [7 7];

FigHandle = f;
figure_width = 8;
figure_height = 6.5;
font_size = 8;
line_width = 1;
export_ppi = 600;
print_png = true;
print_pdf = true;
savename = 'Results/Base/Learning_w_shortage_double';
printsetup(FigHandle, figure_width, figure_height, font_size, line_width, export_ppi, print_png, print_pdf, savename)