%% Run: finite horizion model 

tic

%% Parameters

% Run paramters
runParam = struct;
runParam.runSDP = true;
runParam.simulateOn = true;
runParam.calculateTgw =true;
runParam.saveOn = true; 
runParam.simNum = 1000;
runParam.simpleVersion = false;
runParam.flexOn = true;
runParam.capacityDelay = true;
runParam.solveNoLearning = true;
runParam.adjustOutput = true;
runParam.N = 30;

plotParam = struct;
plotParam.plotsOn = false;
plotParam.policyPlotsOn = true;
plotParam.simPlotsOn = true; 
plotParam.plotInitialWaterBalance = false; 
plotParam.plotHeatMaps = true;
plotParam.plotinfoOverTime = false;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 1;    % $/m^2
% costParam.expansion_cost.capex.large = 258658804 * 2 * .9; % $
% costParam.expansion_cost.capex.small = costParam.expansion_cost.capex.large /3 * 1.15;
costParam.marginal_cost = .48;
costParam.discount_rate = 0.00;

% Population parameters
popParam = struct;
popParam.pop_initial = 6;   % in millions 
popParam.pop_initial = 6.2;
popParam.growth.medium = 0.02;
popParam.growth.high = 0.025;
popParam.growth.low = 0.015;
popParam.growth.none = 0.0;
popParam.growthScenario = 'none';

% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 100;
gwParam.depthLimit = 100;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 17182;
gwParam.wellIndex = 108; % 68 is RR1, 108 is Shemesy, 93 is royal garage
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = true;
gwParam.infoScenario = 'full_range';
gwParam.TgwLoadName = 'T_gw';
gwParam.likelihoodfct = 'normal';
gwParam.llhstddev = 10;

% Water infrastructure paramters
water = struct;
water.desal_capacity_initial = 1.3E6 * 365; % m^3/y
water.desal_capacity_expansion.large = 0.51E6 * 365;
water.desal_capacity_expansion.large = gwParam.pumpingRate * 3;
water.desal_capacity_expansion.small = 0.51E6/3 * 365;
water.desal_capacity_expansion.small = gwParam.pumpingRate;
water.demandFraction = 1;
water.demandPerCapita = 300:-2:300-2*(runParam.N-1);
water.demandPerCapita = 300*ones(1,runParam.N); 
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

% sensInput{1} = deal({'gwParam', 'depthLimit', { 50, 100, 150}});
% sensInput{2} = deal({'gwParam', 'wellIndex', { 108, 93, 68}});
% sensInput{3} = deal({'costParam', 'shortage_cost', { 1, 0.5, 5, 10}});
% sensInput{4} = deal({'costParam', 'discount_rate', { 0, 0.03, 0.05, 0.07}});
% sensInput{5} = deal({'gwParam', 'infoScenario', {'full_range', '10%_cutoff', '20%_cutoff'}});
% sensInput{6} = deal({'gwParam', 'llhstddev', {10, 5}});

sensParams = {...
%     'depthLimit', ...
%     'wellIndex'...
%     'shortage_cost', ...
'discount_rate', ...
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

        [ V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, lowestCostAction, s_gw,...
            s_expand, exp_vectors, K_samples, S_samples, sampleProb] = ...
            sdp_gw( runParam, costParam, popParam, gwParam, water, datetime );
        
        
        evalin('base',  strcat(sensInput{i}{2},'_Output{j} = cell(9,1);'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{1} = sensInput{i}{3}{j};'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{2} = V;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{3} = X1;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{4} = X2;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{5} = T_gw_all;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{6} = lowestCostAction;'));
        evalin('base', strcat(sensInput{i}{2},'_Output{j}{7} = K_samples;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{8} = S_samples;'));
        evalin('base',  strcat(sensInput{i}{2},'_Output{j}{9} = sampleProb;'));
        
        if runParam.saveOn
            save(strcat('sens',sensInput{i}{1}, datetime,'_', num2str(jobid)));
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
    for i = 1 %1:length(sensParams)
        for j = 1:length(sensInput{i})
            
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
    for i = 1 %:length(sensParams)
        figure;
        for j = 1:length(sensInput{i})
            evalin('base', strcat('sim = sim_', sensParams{i}, num2str(j),';'));
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
            yLarge = histc(expTimeLarge,[0:32]);
            ySmall = histc(expTimeSmall,[0:32]);
            subplot(1,length(sensInput{i}), j)
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
            if isnumeric(sensInput{i}{3}{j})
                title(strcat(sensParams{i}, ' = ',  num2str(sensInput{i}{3}{j})));  
            else
                title(strcat(sensParams{i}, ' = ',  sensInput{i}{3}{j}));  
            end
        end
        suptitle(strcat('Histograms of expansion time in ', num2str(R), ' simulations'))
    end
    
    
    % Plot first drawdown level when build (when nocapacity)
    for i = 1:length(sensParams)
        figure;
        for j = 1:length(sensInput{i})
            evalin('base', strcat( 'X2 = ', sensParams{i}, '_Output{j}{4};'));
            if strcmp(sensParams{i}, 'depthLimit')
                gwParam.depthLimit = sensInput{i}{3}{j};
            end
            X2nocap = permute(X2(:,1,1:end-1),[1,3,2]);
            indexZeros = ~flipud(X2nocap == 0);
            indexFirstZero = 200 - sum(cumprod(double(indexZeros),1)) + 3;
            indexNan = sum(cumprod(~isnan(X2nocap)));
            indexReplaceNan = indexFirstZero >= indexNan;
            indexFirstZero(indexReplaceNan) = 1;
            subplot(1,length(sensInput{i}), j)
            scatter(1:N, 200-s_gw(indexFirstZero));
            xlabel('Year')
            ylim([0 205])
            ylabel('Head [m]')
            if isnumeric(sensInput{i}{3}{j})
                title(strcat(sensParams{i}, ' = ',  num2str(sensInput{i}{3}{j})));  
            else
                title(strcat(sensParams{i}, ' = ',  sensInput{i}{3}{j}));  
            end
            hold on
            line([0 30], [200 200], 'Color', 'k')
            line([0 30], [200-gwParam.depthLimit 200-gwParam.depthLimit], 'Color', 'r', 'LineStyle','--')

        end
        suptitle('Optimal expansion policy: Drawdown Threshold for Exapsnion over Time')
    end
    
    
    % Plot total shortage vs total cost
    for i = 1:length(sensParams)
        figure;
        for j = 1:length(sensInput{i})
            evalin('base', strcat('sim = sim_', sensParams{i}, num2str(j),';'));
            totalCost = sum(sim.costOverTime,2);
            totalShortage = sum(sim.shortageOverTime,2);
            subplot(1,length(sensInput{i}), j)
            scatter(totalShortage/gwParam.pumpingRate,totalCost/1E9)
            ylim([0 4])
            xlim([0 8])
            if isnumeric(sensInput{i}{3}{j})
                title(strcat('Avg Cost:', num2str(sim.averageTotalCost, '%.1E'), ',  ', sensParams{i}, ' = ',  num2str(sensInput{i}{3}{j})));  
            else
                title(strcat(sensParams{i}, ' = ',  sensInput{i}{3}{j}));  
            end
        end
        suptitle(strcat('Total shortage vs total cost'))
    end
    
end


%%
toc