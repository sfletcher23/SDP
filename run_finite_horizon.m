%% Run: finite horizion model 

tic

%% Plot parameters


% Run paramters
runSDP = false;
adjustOutput = true;
saveOn = false; % Save output if true
policyPlotsOn = true;
simulateOn = true;
simPlotsOn = true; % Plot results if true
plotInitialWaterBalance = false;
plotHeatMaps = false;
plotSamples = false;
calculateTgw = true;
simpleVersion = false;
infoOverTime = false;
flexOn = true;
capacityDelay = true;

%% Parameters
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
jobid = [];
% turn off plotting if running on cluster
if ~isempty(getenv('SLURM_JOB_ID'))
    disp('job id test working')
    policyPlotsOn = false;
    simPlotsOn = false; % Plot results if true
    plotInitialWaterBalance = false;
    plotHeatMaps = false;
    plotSamples = false;
    jobid = getenv('SLURM_JOB_ID');
end


% Time period
N = 30;

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
gwParam.sampleSize = 1000;
gwParam.depthLimit = 100;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 17182;
gwParam.wellIndex = 108; % 68 is RR1, 108 is Shemesy, 93 is royal garage
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = true;

% Water infrastructure paramters
water = struct;
water.desal_capacity_initial = 1.3E6 * 365; % m^3/y
water.desal_capacity_expansion.large = 0.51E6 * 365;
water.desal_capacity_expansion.large = gwParam.pumpingRate * 3;
water.desal_capacity_expansion.small = 0.51E6/3 * 365;
water.desal_capacity_expansion.small = gwParam.pumpingRate;
water.demandFraction = 1;
water.demandPerCapita = 300:-2:300-2*(N-1);
water.demandPerCapita = 300*ones(1,N); 
if flexOn
    water.desal_capacity_expansion.small = gwParam.pumpingRate/3;
    water.desal_capacity_expansion.large = gwParam.pumpingRate;
else
    water.desal_capacity_expansion.small = gwParam.pumpingRate;
    water.desal_capacity_expansion.large = gwParam.pumpingRate * 3;
end


% Information scenarios
infoScenario = 'full_range';



%% Define population growth and demand
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
if plotInitialWaterBalance
    
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

%% State and Action Definitions for Groundwater 

% Generate state space for groundwater head and demand range
[s_gw, gw_M] = gen_water_growth_states(gwParam);
s_gw = [-1 s_gw]; % This is absorbing state where can't pump anymore
gw_M = gw_M + 1;

% Actions: Stop pumping groundwater (0), continue pumping (1)
a_gw_available = [0 1];
a_gw_unavailable = [0];
  

%% Desalination Expansions: State Definitions and Actions

% a2 desal actions: 0 no expand, 1 expand small, 2 expand large
a_expand_available = [0 1 2];
a_expand_unavailable = 0;

% State definition: volume of additional capacity
maxNumSmallExp = 3;
maxNumLargeExp = 1;

% Check that large capacity is a multiple of small capacity
if mod(water.desal_capacity_expansion.large , water.desal_capacity_expansion.small) ~= 0
    error('Large capacity is not a multiple of small capacity')
end

% Get max capacity, state space between 0 and max cap in steps of small capacity
maxExpCap = water.desal_capacity_expansion.large;
s_expand = 0:water.desal_capacity_expansion.small:maxExpCap;
exp_M = length(s_expand); % Desalination expanded = 2

% Add capacity delay to state space
if capacityDelay
   s_exp_on = s_expand;
   s_exp_delay1 = s_expand;
   s_exp_delay2 = [s_expand(1) s_expand(2)];
   % Index for feasible expansion state combinations - max total 3v across substates
    s_expand = [1:7 9 10 13 17];
    exp_M = length(s_expand);
end


%% Get K and S samples and use to prune state space

if calculateTgw

    [K_samples, S_samples] = gen_param_dist('full_range', gwParam.sampleSize, 1, N);

    % Get min and max hydrograph
    % Get neural net script
    netname = strcat('myNeuralNetworkFunction_', num2str(gwParam.nnNumber));
    netscript = str2func(netname);
    maxK = max(K_samples);
    minK = min(K_samples);
    maxS = max(S_samples);
    minS = min(S_samples);
    time = 1*365:365:N*365;
    x = [ones(1,N) * maxK; ones(1,N) * maxS; time]; 
    y = netscript(x, adjustOutput);
    minDrawdownHydrograph = y(gwParam.wellIndex,:);
    x = [ones(1,N) * minK; ones(1,N) * minS; time]; 
    y = netscript(x, adjustOutput);
    maxDrawdownHydrograph = y(gwParam.wellIndex,:);
    for t = 1:N
        indexValidState = s_gw <= 200 - maxDrawdownHydrograph(t) + 2;
        index_s_gw_time{t} = find(indexValidState);
    end

% Calculate Groudnwater transition matrix when pumping

    % Get transmat vector for gw when pumping for current gw state

    stateInfeasible = true(gw_M, N);
    numRelevantSamples = zeros(gw_M, N);
    indexAbove = cell(gw_M, N);
    indexBelow = cell(gw_M, N);

    T_gw_all = zeros(gw_M, gw_M, N);

    for t =1:N
        parfor index_s1 = 1:gw_M
            s1 = s_gw(index_s1);
            [T_gw_temp, numRel, stateInf, indAbv, indBlw, indRel] = ...
                gw_transrow_nn(gwParam, t, K_samples, S_samples, s1, s_gw, adjustOutput);
            T_gw_all(:,index_s1,t) = T_gw_temp';
            numRelevantSamples(index_s1,t) = numRel;
            stateInfeasible(index_s1,t) = stateInf;
            indexAbove{index_s1, t} = indAbv;
            indexBelow{index_s1, t} = indBlw;
        end
    end    
    
    save(strcat('T_gw_',datetime), 'T_gw_all', 'K_samples', 'S_samples', 'index_s_gw_time')
    
else
    load('T_gw')
end

% Calculate expected total drawdown for each state
cumTgw = zeros(gw_M, N);
for t = linspace(N,1,N)
    for index_s1 = 1:gw_M
        s1 = s_gw(index_s1);
        if t == N 
            cumTgw(index_s1,t) = T_gw_all(1,index_s1,t);
        else
            cumTgw(index_s1,t) = sum(T_gw_all(:,index_s1,t) .* cumTgw(:,t+1));
        end
    end
end

%% Construct simple model version for testing

if simpleVersion
    
s_gw = [-1 0:3]; % This is absorbing state where can't pump anymore
gw_M = length(s_gw);  
N = 4;    
T_gw_all = zeros(gw_M, gw_M, N);
T_gw_all(1,1,:) = 1;
T_gw_all(:,:,1) =  [1    0     0     1/3    1/2;
                    0   1/3    0     0    0;
                    0   1/3    1/3   0    0;
                    0   1/3    1/3   1/3  0;
                    0    0     1/3   1/3  1/2];
T_gw_all(:,:,2) =  [1    0     0     3/5  4/5;
                    0   2/3    0     0    0;
                    0   1/3    1/5   0    0;
                    0    0     3/5   1/5  0;
                    0    0     1/5   1/5  1/5];
T_gw_all(:,:,3) =  [1    0     0     1/5  6/7;
                    0   3/4    0     0    0;
                    0   1/4    2/5   0    0;
                    0    0     2/5   1/5  0;
                    0    0     1/5   3/5  1/7];
T_gw_all(:,:,4) =  [1    0     0     1/7  2/3;
                    0    1     0     0    0;
                    0    0    4/5    0    0;
                    0    0    1/5    1/7  0;
                    0    0     0     5/7  1/3];
                
                
% s_gw = [-1 0 1]; % This is absorbing state where can't pump anymore
% gw_M = length(s_gw);  
% N = 5;    
% T_gw_all = zeros(gw_M, gw_M, N);
% T_gw_all(:,:,1) = [1 0 0;
%                    0 1/2 0
%                    0 1/2 1];
% T_gw_all(:,:,2) = [1 .2 .8;
%                    0 .8 .2
%                    0 0 0];
% T_gw_all(:,:,3) = [1 .1 .9;
%                    0 .9 .1
%                    0 0 0];
% T_gw_all(:,:,4) = [1 .1 .9;
%                    0 .9 .1
%                    0 0 0];
% T_gw_all(:,:,5) = [1 .1 .9;
%                    0 .9 .1
%                    0 0 0];
              

% Check valid p
for t = 1:N
    totprob = sum(T_gw_all(:,:,t),1);
    indexWrong = find(totprob ~= 1);
    sumIndexWrong = sum(indexWrong);
    if indexWrong > 0
        error(strcat('invalid T period ', num2str(t)))
    end
end

for t = 1:N
    index_s_gw_time{t} = 1:gw_M;
end

% Calculate expected total drawdown for each state
cumTgw = zeros(gw_M, N);
for t = linspace(N,1,N)
    for index_s1 = index_s_gw_thisPeriod
        s1 = s_gw(index_s1);
        if t == N 
            cumTgw(index_s1,t) = T_gw_all(1,index_s1,t);
        else
            cumTgw(index_s1,t) = sum(T_gw_all(:,index_s1,t) .* cumTgw(:,t+1));
        end
    end
end



% Cost paramters
costParam = struct;
costParam.shortage_cost = 40;    % $/m^2
% costParam.expansion_cost.capex.large = 258658804 * 2 * .9; % $
% costParam.expansion_cost.capex.small = costParam.expansion_cost.capex.large /3 * 1.15;
costParam.marginal_cost = 2;
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
gwParam.sampleSize = 1000;
gwParam.depthLimit = 0;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 17182;
gwParam.wellIndex = 108; % 68 is RR1, 108 is Shemesy, 93 is royal garage
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = false;

% Water infrastructure paramters
water = struct;
water.desal_capacity_initial = 1.3E6 * 365; % m^3/y
water.desal_capacity_expansion.large = 0.51E6 * 365;
water.desal_capacity_expansion.large = gwParam.pumpingRate * 3;
water.desal_capacity_expansion.small = 0.51E6/3 * 365;
water.desal_capacity_expansion.small = gwParam.pumpingRate;
water.demandFraction = 1;
water.demandPerCapita = 300:-2:300-2*(N-1);
water.demandPerCapita = 300*ones(1,N); 




end

%% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = NaN(gw_M, exp_M, N+1);
X1 = NaN(gw_M, exp_M, N+1);
X2 = NaN(gw_M, exp_M, N+1);

% Terminal period
X1(:,:,N+1) = zeros(gw_M, exp_M, 1);
X2(:,:,N+1) = zeros(gw_M, exp_M, 1);
V(:,:,N+1) = zeros(gw_M, exp_M, 1);

%% Backwards Recursion

% If running on cluster, get number of workers 
if ~isempty(getenv('SLURM_CPUS_PER_TASK'))
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
end

% Loop over all time periods
for t = linspace(N,1,N)
    % Calculate nextV    
    nextV = V(:,:,t+1);
          
    % Loop over all states
    
    % Loop over groundwater state: 1 is depleted, M1 is full
    index_s_gw_thisPeriod = index_s_gw_time{t}; 
    parfor index_s1 = index_s_gw_thisPeriod
        s1 = s_gw(index_s1);
       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for index_s2 = 1:exp_M
            s2 = s_expand(index_s2);
            
            if capacityDelay
                subindex_s2 = linIndex2VecIndex(s2, {s_exp_on', s_exp_delay1', s_exp_delay2'});
            end

            bestV = Inf;
            bestX1= 0;  % Groundwater action and expansion action
            bestX2= 0;
            
            % Update available actions based on whether gw available
            if s1 == 200
                a_gw = a_gw_unavailable;    % unavailble bc depleted
            elseif s1 == -1
                a_gw = a_gw_unavailable;    % unavailble bc turned off
            else
                a_gw = a_gw_available;
            end

            % Update available actions based on whether expansion available
            if capacityDelay
                if subindex_s2(1) == 4 || subindex_s2(2) == 4 || subindex_s2(3) == 2 ... % Max capacity in any of subcategories
                        || (subindex_s2(1) + subindex_s2(2)) >= 5  % Max capacity in online + delay 1
                    a_expand = [0]; % If max capacity online or waiting to come online, can't expand
                elseif (subindex_s2(1) + subindex_s2(2)) == 4 ...   % 2 small units online or in delay 1
                        || (subindex_s2(1) + subindex_s2(2)) == 3    % 1 small units online or in delay 1
                    a_expand = [0 1];
                elseif (subindex_s2(1) + subindex_s2(2) + subindex_s2(3)) == 3
                    a_expand = [0 1 2];
                else
                    error('Some combination was not included!')
                end
                
            else
                switch s2
                    case s_expand(1)
                        a_expand = [0 1 2];
                    case s_expand(2)
                        a_expand = [0 1];
                    case s_expand(3)
                        a_expand = [0 1];
                    case s_expand(4)
                        a_expand = [0];
                end
            end

            num_a_gw = length(a_gw);
            num_a_expand = length(a_expand);

            % Loop over all actions
            % Loop over groundwater pumping action
            for index_a1 = 1:num_a_gw
                a1 = a_gw(index_a1);

                % Loop over expansion action: 1 is do not expand, 2 is expand
                for index_a2 = 1:num_a_expand
                    a2 = a_expand(index_a2);

                    % Calculate demand
%                     demandThisPeriod = demand(water, population(t), t, gwParam);
                    demandThisPeriod = gwParam.pumpingRate;

                    % Calculate cost and shortages this period
                    [ cost, ~,~, ~,~, ~, ~, ~, ~, ~ ] = supplyAndCost( a1, a2, s1, s2, costParam, water, gwParam, t, demandThisPeriod);

                    % Calculate transition matrix
                    
                    % If stop pumping, move to state -1. Otherwise, use
                    % T_gw calculated above. 
                    
                    switch a1
                        case 0
                            T_gw = zeros(1,gw_M);
                            T_gw(1) = 1;
                        case 1
                            T_gw = T_gw_all(:,index_s1,t)';
                    end

                    % Get transmat vector for next expansion state
                    % (deterministic)                  
                    T_expand = zeros(1,exp_M);
                    
                    if capacityDelay
                        T_exp_online_ind = subindex_s2(1);
                        T_exp_delay1_ind = subindex_s2(2);
                        T_exp_delay2_ind = subindex_s2(3);
                        % Move delayed capacity to online
                        if subindex_s2(3) == 2  % Move big plant from delay2 to delay 1
                            T_exp_delay1_ind = 4;
                            T_exp_delay2_ind = 1;
                        elseif subindex_s2(2) > 1
                            T_exp_online_ind = T_exp_online_ind + (subindex_s2(2) - 1);
                            T_exp_delay1_ind = 1;
                        end
                        % Add new capacity to delay
                        if a2 == 1
                            T_exp_delay1_ind = 2;
                        elseif a2 == 2
                            T_exp_delay2_ind = 2;
                        end
                        temp_index = vectorIndex([T_exp_online_ind T_exp_delay1_ind T_exp_delay2_ind], {s_exp_on', s_exp_delay1', s_exp_delay2'});
                        exp_index = find(s_expand == temp_index);
                        T_expand(exp_index) = 1;
                    else
                        if a2 == 0
                            T_expand(index_s2) = 1; % Stay in current state
                        elseif a2 == 1
                            T_expand(index_s2 + 1) = 1; % Move up one state
                        elseif a2 == 2
                            T_expand(index_s2 + 3) = 1; % Move up three states
                        end
                    end
                    
                    % Calculate full transition matrix
                    % T gives probability of next state given
                    % current state and actions

                    TRows = cell(2,1);
                    TRows{1} = T_gw;
                    TRows{2} = T_expand;
                    [ T ] = transrow2mat( TRows );

                     % Calculate expected future cost
                     %nextV = V(:,:,:,:,t+1);
                    indexNonZeroT = find(T > 0);
                    expV = sum(T(indexNonZeroT) .* nextV(indexNonZeroT));
                    for i = 2:4
                        expV = sum(expV);
                    end
                    
                    stateMsg = strcat('t=', num2str(t), ', s1=', num2str(s1), ', a1=', num2str(a1), ', s2=', num2str(s2), ', a2=', num2str(a2))
                    disp(stateMsg)
                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX1 = a1;
                        bestX2 = a2;
                    end
                end
            end
            
%             if saveOn
%                 save(strcat(datetime,'_', num2str(jobid)));
%             end

            % Check that bestV is not Inf
            if bestV == Inf
                error('BestV is Inf, did not pick an action')
            end

            % Save best value and action for current state
            V(index_s1, index_s2, t) = bestV;
            X1(index_s1, index_s2, t) = bestX1;
            X2(index_s1, index_s2, t) = bestX2;

        end
    end
end

if saveOn
    save(strcat(datetime,'_', num2str(jobid)));
end
%% Solve for optimal policies when all decisions made in 1st stage
if false
    
    % Get water demand
    waterDemand = waterDemand_low;
    
    % Get hydrograph for each sample;
    netname = strcat('myNeuralNetworkFunction_', num2str(gwParam.nnNumber));
    netscript = str2func(netname);
    headSample = zeros(length(K_samples), N);
    for i = 1:length(K_samples)
        x = [repmat(K_samples(i),[1,N]); repmat(S_samples(i),[1,N]); [365:365:365*(N)]];
        tempHead = netscript(x, adjustOutput);
        headSample(i,:) = tempHead(gwParam.wellIndex,:);
    end
    
    % Get gw supply and pumping cost for each sample
    gwSupplySample = zeros(length(K_samples), N);
    
end

%% Visualize results: plot optimal policies

if policyPlotsOn
    gw_step = s_gw(3) - s_gw(2);
    exp_step = s_expand(2) - s_expand(1);
    blues = colormap(cbrewer('seq', 'Blues', 6));
    oranges = colormap(cbrewer('seq', 'Oranges', 6));
    color = {blues(2,:), oranges(2,:), blues(4,:), oranges(4,:), blues(6,:), oranges(6,:), [0 0 0]};
    fig = figure;
    times = [22 23 24 25 26 27 28 29 30];
     times = [10 11 12 13 14 15 16];
%     times = [1 2 3 4 5 6 7];
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
                y = [s_expand(j)-(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)-(exp_step/2)];
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
        ax.YTick = s_expand;
        ax.YTickLabel = [];
        xlim([s_gw(1)-gw_step/2 s_gw(end)+gw_step/2])
        ylim([s_expand(1)-exp_step/2 s_expand(end)+exp_step/2])
        if t == 6
            ylabel('Added capacity state [MCM/y]')
            xlabel('Drawdown')
            ax.YTickLabel = string(round(s_expand/1E6));
            ax.XTickLabel = 0:5:s_gw(end);
        end
        title(strcat('Time step: ', num2str(times(t))))
        ax.XTickLabelRotation = 90;
    end
end



%% Plot heat maps
if plotHeatMaps 
    hm1 = HeatMap(flipud(double(~stateInfeasible)), 'Title', 'Infeasible States (in black) from nn samples')
    addXLabel(hm1, 'time')
    addYLabel(hm1, 'drawdown')
    hm2 = HeatMap(flipud(numRelevantSamples), 'Title', 'Number of Relevant Samples: Red is high')
    addXLabel(hm2, 'time')
    addYLabel(hm2, 'drawdown')
    temp = permute(isnan(X1(:,1,:)),[1 3 2]);
    hm3 = HeatMap(flipud(double(~temp)), 'Title', 'Infeasible states (in black) from pruned tree')
end


%% Simulate performance
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

if simulateOn

R = 10;

% Initialize vector tracking state, actions, water balance, costs over time 
state_gw = zeros(R,N);
state_expand = zeros(R,N);
action_gw = zeros(R,N);
action_expand = zeros(R,N);
costOverTime = zeros(R,N);
shortageCostOverTime = zeros(R,N);
expansionCostOverTime = zeros(R,N);
pumpingCostOverTime = zeros(R,N);
shortageOverTime = zeros(R,N);
capacityOverTime = zeros(R,N);
minjurSupplyOverTime = zeros(R,N);
othergwSupplyOverTime = zeros(R,N);
demandOverTime = zeros(R,N);
expSupplyOverTime = zeros(R,N);
margDesalCostOverTime = zeros(R,N);
T_gw_time = zeros(gw_M,N,R);
sampleIndexOverTime = zeros(gwParam.sampleSize,N,R);
failureProbOverTime = zeros(R,N);

% Initial state
s_gw_initial = 0;
s_expand_initial = 0;

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;

[K_samples, S_samples] = gen_param_dist(infoScenario, gwParam.sampleSize, 1, N);

for i = 1:R
    
    state_gw_now = zeros(1,N);
    state_expand_now = zeros(1,N);
    action_gw_now = zeros(1,N);
    action_expand_now = zeros(1,N);
    costOverTime_now = zeros(1,N);
    shortageCostOverTime_now = zeros(1,N);
    expansionCostOverTime_now = zeros(1,N);
    pumpingCostOverTime_now = zeros(1,N);
    shortageOverTime_now = zeros(1,N);
    capacityOverTime_now = zeros(1,N);
    minjurSupplyOverTime_now = zeros(1,N);
    othergwSupplyOverTime_now = zeros(1,N);
    demandOverTime_now = zeros(1,N);
    expSupplyOverTime_now = zeros(1,N);
    margDesalCostOverTime_now = zeros(1,N);
    T_gw_time_now = zeros(gw_M,N);
    sampleIndexOverTime_now = zeros(gwParam.sampleSize,N);
    failureProbOverTime_now = zeros(1,N);
    
    for t = 1:N

        % Caculate state indexes
        index_state_gw = find(state_gw_now(t) == s_gw);
        index_state_expand = find(state_expand_now(t) == s_expand);
        
        % Lookup failure prob if keep pumpin
        failureProbOverTime_now(t) = cumTgw(index_state_gw ,t);
        
        % Lookup optimal policy for current state
        action_gw_now(t) = X1(index_state_gw, index_state_expand, t);
        action_expand_now(t) = X2(index_state_gw, index_state_expand, t);
%         if state_gw_now(t) < gwParam.depthLimit
%             action_gw_now(t) = 1;
%         end
%         if t < 5
%             action_expand_now(t) = 0;
%         elseif t == 5
%             action_expand_now(t) = 1;
%         end
%          action_expand_now(t) = 0;    
%            if t ==1 
%                 action_expand_now(t) = 1;
%            end
        % Calculate demand, shortage, and cost for current t
        demandOverTime_now(t) = demand( water, population(t), t, gwParam);
        [ costOverTime_now(t), shortageCostOverTime_now(t), expansionCostOverTime_now(t), pumpingCostOverTime_now(t), margDesalCostOverTime_now(t), ...
            shortageOverTime_now(t), capacityOverTime_now(t),minjurSupplyOverTime_now(t), expSupplyOverTime_now(t), othergwSupplyOverTime_now ] ...
            = supplyAndCost( action_gw_now(t),  action_expand_now(t), state_gw_now(t), state_expand_now(t), costParam, water, gwParam, t,  demandOverTime_now(t));
        
        % Get transisition mat to next state give current state and actions

            % Get transmat vector to next GW state 
            if action_gw_now(t) == 0
                T_current_gw = zeros(1, length(s_gw));
                T_current_gw(1) = 1;
                index = ones(gwParam.sampleSize,1)*-99;
            else
                T_current_gw = T_gw_all(:,index_state_gw,t)';
            end
%             T_gw_time_now(:,t) = T_current_gw;
%             sampleIndexOverTime_now(:,t) = index;

            % Get transmat vector for next expansion state (deterministic)
            T_current_expand = zeros(1,exp_M);
            if action_expand_now(t) == 0
                T_current_expand(index_state_expand) = 1; % Stay in current state
            elseif action_expand_now(t) == 1
                T_current_expand(index_state_expand + 1) = 1; % Move up one state
            elseif action_expand_now(t) == 2
                T_current_expand(index_state_expand + 3) = 1; % Move up three states
            end

            % Get Transition Matrix from rows
            TRows_current = cell(2,1);
            TRows_current{1} = T_current_gw;
            TRows_current{2} = T_current_expand;
            [ T_current ] = transrow2mat( TRows_current );

        % Simulate next state
        if t < N
            T_current_1D = reshape(T_current,[1 numel(T_current)]);
            T_current_1D_cumsum = cumsum(T_current_1D);
            p = rand();
            index = find(p < T_current_1D_cumsum,1);
            [ind_s1, ind_s2, ind_s3, ind_s4] = ind2sub(size(T_current),index);
                % Test sample
                margin = 1e-10;
                if (T_current(ind_s1, ind_s2, ind_s3, ind_s4) < margin)
                    error('Invalid sample from T_current')
                end

            state_gw_now(t+1) = s_gw(ind_s1); 
            state_expand_now(t+1) = s_expand(ind_s2);
                % Test next state
                test_gw = T_current_gw(ind_s1) >= 0;
                test_expand = T_current_expand(ind_s2) >= 0;
                if ~test_gw
                    error('Invalid gw state tranisition')
                end
                if ~test_expand
                    error('Invalid expand state tranisition')
                end
        end
    end
    state_gw(i,:) = state_gw_now;
    state_expand(i,:) = state_expand_now;
    action_gw(i,:) = action_gw_now;
    action_expand(i,:) = action_expand_now;
    costOverTime(i,:) = costOverTime_now;
    shortageCostOverTime(i,:) = shortageCostOverTime_now;
    expansionCostOverTime(i,:) = expansionCostOverTime_now;
    pumpingCostOverTime(i,:) =pumpingCostOverTime_now;
    shortageOverTime(i,:) = shortageOverTime_now;
    capacityOverTime(i,:) = capacityOverTime_now;
    minjurSupplyOverTime(i,:) = minjurSupplyOverTime_now;
    othergwSupplyOverTime(i,:) = othergwSupplyOverTime_now;
    demandOverTime(i,:) = demandOverTime_now;
    expSupplyOverTime(i,:) = expSupplyOverTime_now;
    margDesalCostOverTime(i,:) = margDesalCostOverTime_now;
    T_gw_time(:,:,i) = T_gw_time_now;
    sampleIndexOverTime(:,:,i) = sampleIndexOverTime_now;
    failureProbOverTime(i,:) = failureProbOverTime_now;
end

T_gw_save = T_gw_time_now;

end

if simPlotsOn

if R == 1

    % Plot state evolution w/ actions
    figure;
    yyaxis left
    plot(1:N, state_gw)
    hold on
    yyaxis right
    plot(1:N, action_gw)
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
    plot(1:N, failureProbOverTime)


    subplot(1,3,2)
    plot(1:N,costOverTime/1E6);
    h = gca;
    % ylim([0 700])
    hold on
    bar(1:N, [shortageCostOverTime./1E6; expansionCostOverTime./1E6; pumpingCostOverTime./1E6; margDesalCostOverTime./1E6]', 'stacked');
    legend('Total cost', 'Shortage cost', 'Expansion Cost', 'Pumping Cost', 'Desal costs')
    title(strcat('Total cost [M$]: ', num2str(sum(costOverTime)/1E6, '%.3E')))
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
    plot(1:N,shortageOverTime/1E6)
    hold on
    plot(1:N,demandOverTime/1E6)
    plot(1:N,capacityOverTime/1E6)
    bar(1:N, [minjurSupplyOverTime/1E6; expSupplyOverTime/1E6]', 'stacked');
    legend('shortage', 'demand', 'capacity', 'minjur supply', 'exp supply')
    legend('Location', 'southwest')
    ylabel('MCM/y');
    
else
    % Plot hydrographs
    figure;
    indexLimit = find(state_gw == -1 | state_gw >= 100);
    hydrograph = 200 - state_gw;
    hydrograph(indexLimit) = gwParam.depthLimit;
    plot(1:N, hydrograph)
    xlabel('Time [years]')
    ylabel('Drawdown [m]')
    title('Simulated Drawdown')
    
    % Plot expansion time distribution
    [r,c] = find(expansionCostOverTime);
    % for every row, take the minimum column index and put NaN if none is found
    expTime = accumarray(r,c,[size(expansionCostOverTime,1),1],@min,32);
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
    totalCost = sum(costOverTime,2);
    totalShortage = sum(shortageOverTime,2);
    figure
    scatter(totalShortage,totalCost)
    
    
    % Bagplot!! Eventually, have different bags for learning vs no learning
    % and flexible vs no flexible
    
end

end

%% Show updated predictions over time
if infoOverTime
    
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
                tempHead = netscript(x, adjustOutput);
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

%% Save results

if saveOn
    save(strcat(datetime,'_', num2str(jobid)));
end


toc
