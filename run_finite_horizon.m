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
costParam.discount_rate = 0.05;


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

%% Run SDP  

% Run SDP model. Four subparts can run: calculating T_gw, running simple
% version of SDP, running full version of SDP, choosing best option when
% restricted to 1st period only. 

if runParam.runSDPfunction

    % If running on cluster, get number of workers 
    if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
        parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
    end
    
        [ V, X1, X2, T_gw_all, cumTgw, s_gw, s_expand, exp_vectors, lowestCost, lowestCostAction ] = ...
        sdp_gw( runParam, costParam, gwParam, water, datetime );

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
        lowestCostAction = 2;
        [ simnolearn_build ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
        lowestCostAction = 0;
        [ simnolearn_nobuild ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
    end
    
    if runParam.saveOn
    save(strcat(datetime,'_', num2str(jobid)));
    end
    
end


%% Make plots

if plotParam.plotsOn
	plots_sdp_gw( V, X1, X2, T_gw_all, cumTgw, lowestCost, lowestCostAction, sim, simnolearn_build, simnolearn_nobuild, plotParam, s_gw, s_expand, exp_vectors, runParam, gwParam, costParam, water);
    
end


%%
toc
