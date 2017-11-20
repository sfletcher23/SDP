%% Run: finite horizion model 

tic

%% Parameters


% Run paramters
runParam = struct;
runParam.runSDP = true;
runParam.simulateOn = false;
runParam.calculateTgw = true;
runParam.saveOn = true; 
runParam.simNum = 100;
runParam.simpleVersion = false;
runParam.flexOn = true;
runParam.capacityDelay = true;
runParam.solveNoLearning = true;
runParam.adjustOutput = true;
runParam.runSDPfunction = true;
runParam.oldCost = true;
runParam.N = 30;

plotParam = struct;
plotParam.plotsOn = false;
plotParam.policyPlotsOn = true;
plotParam.simPlotsOn = true; 
plotParam.plotInitialWaterBalance = false; 
plotParam.plotHeatMaps = false;
plotParam.plotinfoOverTime = false;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 10    ;    % $/m^2
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
gwParam.sampleSize = 2000;
gwParam.depthLimit = 100;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 54215;
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = true;
gwParam.infoScenario = 'full_range';
gwParam.TgwLoadName = 'T_gw';
gwParam.likelihoodfct = 'normal';
gwParam.llhstddev = 10;
gwParam.startingHead = 337.143;
gwParam.nstp = 100;

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

%% Run SDP  

% Run SDP model. Four subparts can run: calculating T_gw, running simple
% version of SDP, running full version of SDP, choosing best option when
% restricted to 1st period only. 

if runParam.runSDPfunction

    % If running on cluster, get number of workers 
    if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
        parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
    end
    
    [ V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, lowestCostAction, ...
        s_gw, s_expand, exp_vectors, K_samples, S_samples, sampleProb ] = ...
        sdp_gw( runParam, costParam, popParam, gwParam, water, datetime );

    if runParam.saveOn
        save(strcat(datetime,'_', num2str(jobid)));
    end
    
end
    
%% Run forward simulation 
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

if runParam.simulateOn
% 
% poolobj = gcp;
% addAttachedFiles(poolobj,{'supplyAndCost.m'})

    useNoInfoPolicy = false;
    [ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

    if runParam.solveNoLearning
        useNoInfoPolicy = true;
        [ simnolearn ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
    end
    
    if runParam.saveOn
    save(strcat(datetime,'_', num2str(jobid)));
    end
    
end


%% Make plots

if plotParam.plotsOn
	plots_sdp_gw(  V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, ...
        lowestCostAction, sim, plotParam, s_gw, s_expand, exp_vectors, runParam, gwParam, costParam, water, sampleProb, K_samples, S_samples);
    
    if runParam.solveNoLearning
        plots_sdp_gw(  V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, ...
        lowestCostAction, simnolearn, plotParam, s_gw, s_expand, exp_vectors, runParam, gwParam, costParam, water, sampleProb, K_samples, S_samples);
    
    % Combined plot
    figure;
    totalCost = sum(sim.costOverTime,2);
    totalShortage = sum(sim.shortageOverTime,2);
    scatter(totalShortage/1E6,totalCost/1E9, 50)
    totalCost = sum(simnolearn.costOverTime,2);
    totalShortage = sum(simnolearn.shortageOverTime,2);
    hold on
    scatter(totalShortage/1E6,totalCost/1E9)
    title('Total Shortage vs. Total Cost')
    legend({strcat('Learning over time policy: Average Cost ', num2str(sim.averageTotalCost, '%.3E')), ...
        strcat('Fixed policy: Average Cost ', num2str(simnolearn.averageTotalCost, '%.3E'))},'FontSize', 12)
    ylim([0 5])
    xlabel('Total Shortage over 30 years [MCM]')
    ylabel('Total Costs (including shortages) over 30 years [Billion USD]')
    end

end


%%
toc
