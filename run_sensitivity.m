%% Run: finite horizion model 

tic

%% Parameters


% Run paramters
runParam = struct;
runParam.runSDP = true;
runParam.simulateOn = true;
runParam.calculateTgw = false;
runParam.saveOn = true; 
runParam.simNum = 1000;
runParam.simpleVersion = false;
runParam.flexOn = true;
runParam.capacityDelay = true;
runParam.solveNoLearning = true;
runParam.adjustOutput = true; 
runParam.runSDPfunction = true;
runParam.oldCost = false;
runParam.percentile = 0;
runParam.N = 30;

plotParam = struct;
plotParam.plotsOn = true;
plotParam.policyPlotsOn = false;
plotParam.simPlotsOn = true; 
plotParam.plotInitialWaterBalance = false; 
plotParam.plotHeatMaps = false;
plotParam.plotinfoOverTime = true;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 20;    % $/m^2
costParam.marginal_cost = .48;
costParam.discount_rate = 0.0;


% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 2000;
gwParam.depthLimit = 50;
gwParam.pumpingRate = 298254 * 365;  % m^3/y
gwParam.nnNumber = 54212;
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = true;
gwParam.infoScenario = 'full_range';
gwParam.TgwLoadName = 'T_gw_Dec5';
gwParam.startingHead = 337.143;
gwParam.nstp = 100;

% Water infrastructure paramters
water = struct;
if runParam.flexOn
    water.desal_capacity_expansion.small = gwParam.pumpingRate/3;
    water.desal_capacity_expansion.large = gwParam.pumpingRate;
else
    water.desal_capacity_expansion.small = gwParam.pumpingRate;
    water.desal_capacity_expansion.large = gwParam.pumpingRate * 3;
end

% Set up saving
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
jobid = [];

% Turn off plotting if running on cluster
if ~isempty(getenv('SLURM_JOB_ID'))
    disp('job id test working')
    plotParam.PlotsOn = false;
    plotParam.policyPlotsOn = false;
    plotParam.simPlotsOn = false; % Plot results if true
    plotParam.plotInitialWaterBalance = false;
    plotParam.plotHeatMaps = false;
    plotParam.plotSamples = false;
    jobid = getenv('SLURM_JOB_ID');
end

%% Sensitivity inputs

% sensInput = cell(6,1);

sensInput{1} = deal({'costParam', 'discount_rate', { 0, 0.03, 0.05, 0.07}});
sensInput{2} = deal({'costParam', 'shortage_cost', { 5, 10, 20, 50}});

% sensInput{1} = deal({'gwParam', 'depthLimit', { 50, 100, 150}});
% sensInput{4} = deal({'costParam', 'discount_rate', { 0, 0.03, 0.05, 0.07}});
% sensInput{5} = deal({'gwParam', 'infoScenario', {'full_range', '10%_cutoff', '20%_cutoff'}});
% sensInput{6} = deal({'gwParam', 'llhstddev', {10, 5}});

sensParams = {...
%     'depthLimit', ...
%     'wellIndex'...
'discount_rate' ,...
'shortage_cost' ...
% 'infoScenario', 'llhstddev'...
};



%% Run SDP  

% Run SDP model. Four subparts can run: calculating T_gw, running simple
% version of SDP, running full version of SDP, choosing best option when
% restricted to 1st period only. 

% If running on cluster, get number of workers 
if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
end

for i = 1:length(sensInput)
    
    % Initialize output
%     sens.(sensInput{i}{2}) = cell(length(sensInput{i}{3}),1);
    evalin('base', strcat(sensInput{i}{2},'_Output', '=', 'cell(length(sensInput{i}{3}),1);'));
    
    for j = 1:length(sensInput{i}{3})
        
        % Change input parameter for sensitivity
        if isnumeric(sensInput{i}{3}{j})
            evalin('base', strcat(sensInput{i}{1},'.',sensInput{i}{2}, '=',num2str(sensInput{i}{3}{j}),';'));
        else
            evalin('base', strcat(sensInput{i}{1},'.',sensInput{i}{2}, '='' ', sensInput{i}{3}{j}, '''',';'));
        end

        [ V, X1, X2, T_gw_all, cumTgw, s_gw, s_expand, exp_vectors, lowestCost, lowestCostAction] = ...
            sdp_gw( runParam, costParam, gwParam, water, datetime );
        
        
        evalin('base',  strcat(sensInput{i}{2},'_Output{j} = cell(9,1);'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{1} = sensInput{i}{3}{j};'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{2} = V;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{3} = X1;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{4} = X2;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{5} = T_gw_all;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{6} = lowestCostAction;'));
        
        if runParam.saveOn
            %save(strcat('sens',sensInput{i}{1}, datetime,'_', num2str(jobid)));
        end
        
        % Change parameter back to base case for next round
        if isnumeric(sensInput{i}{3}{j})
            evalin('base', strcat(sensInput{i}{1},'.',sensInput{i}{2}, '=',num2str(sensInput{i}{3}{1})));
        else
            evalin('base', strcat(sensInput{i}{1},'.',sensInput{i}{2}, '='' ', sensInput{i}{3}{1}, ''''));
        end
    end
end
%% Run forward simulation 
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

if runParam.simulateOn
    for i = 1:length(sensParams)
        for j = 1:length(sensInput{i}{3})
            
            evalin('base', strcat( 'V = ', sensParams{i}, '_Output{j}{2};'));
            evalin('base', strcat( 'X1 = ', sensParams{i}, '_Output{j}{3};'));
            evalin('base', strcat( 'X2 = ', sensParams{i}, '_Output{j}{4};'));
            evalin('base', strcat( 'T_gw_all = ', sensParams{i}, '_Output{j}{5};'));
            evalin('base', strcat( 'lowestCostAction = ', sensParams{i}, '_Output{j}{6};'));
            
            useNoInfoPolicy = false;

            [ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

            evalin('base', strcat('sim_', sensParams{i}, num2str(j), ' = sim;'));
            
            if runParam.saveOn
            save(strcat(datetime,'_sim_',sensParams{i}, num2str(jobid)));
            end
        end
    end
end


%% Make plots

if plotParam.plotsOn
    R = runParam.simNum;
    N = 30;
    
    % Expansion time sensitivity
    fig = cell(1,4);
    xlabels = cell(1,4);
    ylabels = cell(1,4);
    xticks = cell(1,4);
    leg = cell(1,4);
    axs = cell(1,4);
    for i = 1 :length(sensParams)
        fig{i} = figure;
        xticktext = cell(1,length(sensInput{i}{3}));
        for j = 1:length(sensInput{i}{3})
            evalin('base', strcat('sim = sim_', sensParams{i}, num2str(j),';'));
            % Plot expansion time distribution
            [~,~, largeCost,~,~,~,~,~,~] = supplyAndCost( 0, 2, 0, 0, costParam, water, gwParam, 1, gwParam.pumpingRate, runParam.capacityDelay, exp_vectors, false);
            indexLarge = sim.state_expand == 17;
            expLargeOverTime = zeros(size(sim.state_expand));
            expLargeOverTime(indexLarge) = 1;
            [rLarge,cLarge] = find(expLargeOverTime);
            % for every row, take the minimum column index and put NaN if none is found
            expTimeLarge = accumarray(rLarge,cLarge,[size(expLargeOverTime,1),1],@min,32);
            countNever = sum((expTimeLarge == 32 ));
            yLarge = histc(expTimeLarge, [0:32]);
            dec1exp = sum(yLarge(1:10));
            dec2exp = sum(yLarge(11:20));
            dec3exp = sum(yLarge(21:30));
            nevexp = yLarge(end);
            bardata = zeros(4);
            bardata(j,:) = [dec1exp dec2exp dec3exp nevexp];
            bar(1:4, bardata, 'stacked')
            hold on 
            xticktext{j} = sensInput{i}{3}{j};
        end
        axs{i} = gca;
        axs{i}.XTick = 1:4;
        axs{i}.XTickLabel = xticktext;
        xticks{i} = xticktext;
        xlabels{i} = xlabel('Expansion Year');
        ylabels{i} = ylabel('Frequency (in 1000 simulations)');
        title(strcat('Sensitivity of expansion time to',{' '}, strrep(sensParams{i},'_', ' ')))
        leg{i} = legend('1st Decade', '2nd Decade', '3rd Decade', 'Never');
    end
    
    
    % Plot first drawdown level when build (when nocapacity)
    for i = 1:length(sensParams)
        fig{i+2} = figure;
        legText = cell(1,length(sensInput{i}{3})+2);
        for j = 1:length(sensInput{i}{3})
            evalin('base', strcat( 'X2 = ', sensParams{i}, '_Output{j}{4};'));
            if strcmp(sensParams{i}, 'depthLimit')
                gwParam.depthLimit = sensInput{i}{3}{j};
            end
            
            X2nocap = permute(X2(:,1,1:end-1),[1,3,2]);
            X2nocap(1,:) = 0;
            indexNan = isnan(X2nocap);
            X2nocap(indexNan) = 0;
            indexZeros = ~flipud(X2nocap == 0);
            indexFirstZero =length(s_gw) - sum(cumprod(double(indexZeros),1));
            indexNan = sum(cumprod(~isnan(X2nocap)));
            indexReplaceNan = indexFirstZero >= indexNan;
            indexFirstZero(indexReplaceNan) = 1;
            if j == 1
                line([0 30], [gwParam.startingHead, gwParam.startingHead ], 'Color', 'k')
                line([0 30], [gwParam.startingHead-gwParam.depthLimit gwParam.startingHead-gwParam.depthLimit], 'Color', 'r', 'LineStyle','--')
                legText{1} = 'Starting head';
                legText{2} = 'Depth limit';
            end
            hold on
            plot(1:N, gwParam.startingHead -   s_gw(indexFirstZero), '-o');
            xlabels{i+2} = xlabel('Year');
            ylabels{i+2} = ylabel('Head [m]');
            hold on
            legText{j+2} = num2str(sensInput{i}{3}{j});
            axs{i+2} = gca;
        end
        title(strcat('Sensitivity of drawdown threshold for exapsnion to ',{' '}, strrep(sensParams{i},'_', ' ')))
        leg{i+2} = legend(legText);
    end
    
    
    % Plot total cost cdf
    if false
    for i = 1:length(sensParams)
        figure;
        legtext = cell(1,length(sensInput{i}{3}));
        for j = 1:length(sensInput{i}{3})
            evalin('base', strcat('sim = sim_', sensParams{i}, num2str(j),';'));
            totalCost = sum(sim.costOverTime,2);
            c = cdfplot(totalCost/1E9);
            c.LineWidth = 1.5;
            hold on 
            legtext{j} = num2str(sensInput{i}{3}{j});
        end
        title(strcat('Sensitivity of cost (with damages) to',{' '}, strrep(sensParams{i},'_', ' ')))
        legend(legtext)
    end
    end
    
    
    % Combine 
    f = figure;
    for i=1:4
        sub(i) = subplot(2,2,i);
    end
    for i=1:4
        copyobj([leg{i} axs{i}],f);
        copyobj(allchild(get(fig{i},'CurrentAxes')),sub(i));
        copyobj(xlabels{i},sub(i));
        copyobj(ylabels{i},sub(i));
    end
    for i=1:2
        copyobj(xticks{i},sub(i));
    end
end



%%
toc