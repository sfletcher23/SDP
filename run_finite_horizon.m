%% Run: finite horizion model 

tic

%test comment
%% Parameters

% Run paramters
policyPlotsOn = true;
parforOn = false; % Parallel processing on?
simulateOn = true;
simPlotsOn = true; % Plot results if true
saveOn = false; % Save output if true
plotInitialWaterBalance = true;

% Time period
N = 30;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 50;
costParam.expansion_cost = 10000000; 
costParam.pumping_cost = 1000;
costParam.discount_rate = 0.04;

% Water infrastructure paramters
water = struct;
water.desal_capacity_initial = 1.3E6 * 365; % m^3/y
water.desal_capacity_expansion = 0.5E6 * 365;
water.demandFraction = 1;
water.demandPerCapita = 300:-2:300-2*(N-1);

% Population parameters
popParam = struct;
popParam.pop_initial = 6;   % in millions 
popParam.growth.medium = 0.03;
popParam.growth.high = 0.04;
popParam.growth.low = 0.02;
popParam.growthScenario = 'medium';

% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 10000;
gwParam.depthLimit = 200;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 17182;
gwParam.wellIndex = 55;


% Information scenarios
infoScenario = 'high_narrow';



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
    population_low(1) = popParam.pop_initial;
    population_medium(1) = popParam.pop_initial;
    population_high(1) = popParam.pop_initial;
    for t = 2:N
        growthRate = popParam.growth.low;
        population_low(t) = population_low(t-1) * (1 + growthRate);
        growthRate = popParam.growth.medium;
        population_medium(t) = population_medium(t-1) * (1 + growthRate);
        growthRate = popParam.growth.high;
        population_high(t) = population_high(t-1) * (1 + growthRate);
    end
    
    gw_Minjur = ones(1,N) * gwParam.pumpingRate;
    gw_other = ones(1,N) * gwParam.otherPumpingRate;
    desal = ones(1,N) * water.desal_capacity_initial;
    desal_exp = ones(1,N) * water.desal_capacity_expansion; 
    waterDemand_low = demand(water, population_low, 1:N);
    waterDemand_medium = demand(water, population_medium, 1:N);
    waterDemand_high = demand(water, population_high, 1:N);
    figure;
    subplot(1,2,1)
    area(1:N, [gw_Minjur; gw_other; desal; desal_exp]' ./ 1E6);
    hold on;
    plot(1:N, waterDemand_low/1E6)
    plot(1:N, waterDemand_medium/1E6)
    plot(1:N, waterDemand_high/1E6)
    legend('Minjur GW', 'Other GW', 'Desal', 'Desal Expansion', 'Demand Low', 'Demand Medium', 'Demand High')
    legend('Location','northwest')
    ylabel('MCM/y')
    xlabel('Year')
    title('Water Balance: With Minjur')
    subplot(1,2,2)
    area(1:N, [gw_other; desal; desal_exp]' ./ 1E6);
    hold on;
    plot(1:N, waterDemand_low/1E6)
    plot(1:N, waterDemand_medium/1E6)
    plot(1:N, waterDemand_high/1E6)
    legend('Other GW', 'Desal', 'Desal Expansion', 'Demand Low', 'Demand Medium', 'Demand High')
    legend('Location','northwest')
    ylabel('MCM/y')
    xlabel('Year')
    title('Water Balance: Without Minjur')
    
end

%% State and Action Definitions for Groundwater 

% Generate state space for groundwater head and demand range
[s_gw, gw_M] = gen_water_growth_states(gwParam);

% Actions: Pump groundwater this period (full demand) or not
    % 1 is pump, 0 is no pump
a_gw_available = [0 1];
a_gw_depleted = [0];


%% Desalination Expansions: State Definitions and Actions

% State definitions
s_expand = 1:2;
exp_M = length(s_expand); % Desalination expanded = 2

% a2 desal actions: 0 no expand, 1 expand
a_expand_available = [0 1];
a_expand_unavailable = [0];


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
if exist(getenv('SLURM_CPUS_PER_TASK'))
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
end

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Calculate nextV    
    nextV = V(:,:,t+1);
    
    % Get K S samples for this period
    [K_samples_thisPeriod, S_samples_thisPeriod] = gen_param_dist(infoScenario, gwParam, t, N);
        
    % Loop over all states
    
    % Loop over groundwater state: 1 is depleted, M1 is full
    parfor index_s1 = 1:gw_M 
        s1 = s_gw(index_s1);
       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for index_s2 = 1:exp_M
            s2 = s_expand(index_s2);

            bestV = Inf;
            bestX = [0 0];  % Groundwater action and expansion action

            % Update available actions based on whether gw depleted
            if s1 == max(s_gw) 
                a_gw = a_gw_depleted;
            else
                a_gw = a_gw_available;
            end

            % Update available actions based on whether expansion available
            if s2 < 2
                a_expand = a_expand_available;
            else
                a_expand = a_expand_unavailable;
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
                    demandThisPeriod = demand(water, population(t), t);

                    % Calculate cost and shortages this period
                    [shortage, ~, ~, gw_supply] =  shortageThisPeriod(a1, s1, s2, water, demandThisPeriod, s_gw, gwParam);
                    cost = costThisPeriod(a1, a2, costParam, shortage, gw_supply, t);

                    % Calculate transition matrix

                    % Get transmat vector for gw based on action, current
                    % gw state
                    T_gw = gw_transrow_nn(gwParam.nnNumber, gwParam.wellIndex, t, K_samples_thisPeriod, S_samples_thisPeriod, s1, s_gw  );

                    % Get transmat vector for next expansion state
                    % (deterministic)
                    if a2 == 1 || s2 == 2   % desal already expanded or will expand
                        T_expand = [0 1];
                    else
                        T_expand = [1 0];
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

                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX = [a1 a2];
                    end

                    bestX1 = a1;
                    bestX2 = a2;
                end
            end

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

%% Visualize results: plot optimal policies

if policyPlotsOn
    gw_step = s_gw(2) - s_gw(1);
    exp_step = s_expand(2) - s_expand(1);
    color = {'b', 'g', 'y', 'r'};
    fig = figure;
    for t = 1:N
        subplot(N,1,t)
        if t == 1
            patch(1,1,color{1}) % Just to make legend work, will be covered up later
            patch(1,1,color{2})
            patch(1,1,color{3})
            patch(1,1,color{4})
            leg = legend('No pump, no expand', 'Pump, no expand', 'No pump, expand', 'Pump, expand');
%                 leg.Location = 'southeastoutside';
        end
        for i = 1:gw_M
            for j = 1:exp_M
                x = [s_gw(i)-(gw_step/2) s_gw(i)-(gw_step/2) s_gw(i)+(gw_step/2) s_gw(i)+(gw_step/2)];
                y = [s_expand(j)-(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)-(exp_step/2)];
                if X1(i,j,1) == 0 && X2(i,j,1) == 0
                    colorThisState = color{1};
                elseif X1(i,j,1) == 1 && X2(i,j,1) == 0
                    colorThisState = color{2};
                elseif X1(i,j,1) == 0 && X2(i,j,1) == 1
                    colorThisState = color{3};
                elseif X1(i,j,1) == 1 && X2(i,j,1) == 1
                     colorThisState = color{4};
                end
                patch(x,y,colorThisState)
                hold on  
            end
        end
    
    ax = gca;
    ax.XTick = 0:5:s_gw(end);
    ax.YTick = s_expand;
    xlim([s_gw(1)-gw_step/2 s_gw(end)+gw_step/2])
    ylim([s_expand(1)-exp_step/2 s_expand(end)+exp_step/2])
    xlabel('Groundwater state')
    ylabel('Expand state')
    ax.YTickLabel = {'Not expanded', 'Expanded'};
    title(strcat('Time step: ', num2str(t)))
    ax.XTickLabelRotation = 90;
    end
end


%% Simulate performance
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

if simulateOn

% Initialize vector tracking state, actions, water balance, costs over time 
state_gw = zeros(1,N);
state_expand = zeros(1,N);
action_gw = zeros(1,N);
action_expand = zeros(1,N);
costOverTime = zeros(1,N);
shortageCostOverTime = zeros(1,N);
expansionCostOverTime = zeros(1,N);
pumpingCostOverTime = zeros(1,N);
shortageOverTime = zeros(1,N);
supplyOverTime = zeros(1,N);
gwSupplyOverTime = zeros(1,N);
demandOverTime = zeros(1,N);

% Initial state
s_gw_initial = s_gw(1);
s_expand_initial = 1;

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;

for t = 1:N
    
    % Caculate state indexes
    index_state_gw = find(state_gw(t) == s_gw);
    index_state_expand = find(state_expand(t) == s_expand);
    
    % Lookup optimal policy for current state
    action_gw(t) = X1(index_state_gw, index_state_expand, t);
    action_expand(t) = X2(index_state_gw, index_state_expand, t);
    
    % Calculate demand, shortage, and cost for current t
    demandOverTime(t) = demand( water, population(t), t);
    [shortageOverTime(t), supplyOverTime(t), ~, gwSupplyOverTime(t)] = shortageThisPeriod(action_gw(t), ...   
        state_gw(t), state_expand(t), water, demandOverTime(t), s_gw, gwParam);
    [costOverTime(t), shortageCostOverTime(t), expansionCostOverTime(t), pumpingCostOverTime(t)]  = ...
        costThisPeriod(action_gw(t), action_expand(t), costParam, shortageOverTime(t),gwSupplyOverTime(t),t);  
    
    % Get transisition mat to next state give current state and actions

        % Get transmat vector to next GW state 
        [K_samples_thisPeriod, S_samples_thisPeriod] = gen_param_dist(infoScenario, gwParam, t, N);
        T_current_gw = gw_transrow_nn(gwParam.nnNumber, gwParam.wellIndex, t, K_samples_thisPeriod, S_samples_thisPeriod, state_gw(t), s_gw );     
 
        % Get transmat vector for next expansion state (deterministic)
        if action_expand(t) == 1 || state_expand(t) == 2   % desal already expanded or will expand
            T_current_expand = [0 1];
        else
            T_current_expand = [1 0];
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
            
        state_gw(t+1) = s_gw(ind_s1); 
        state_expand(t+1) = s_expand(ind_s2);
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

end

%% Plot simulation results

if simPlotsOn

% Plot state evolution w/ actions
figure;
subplot(2,2,1)
yyaxis left
plot(1:N, state_gw')
hold on
yyaxis right
scatter(1:N, action_gw)
xlabel('time')
legend('GW state', 'pumping level')

subplot(2,2,2)
plot(1:N, state_expand')
hold on
scatter(1:N, action_expand')
xlabel('time')
legend('Expansion state', 'Expansion decision')


% Plot system performance
subplot(2,2,3)
plot(1:N,costOverTime);
h = gca;
h.YLim(1) = 0;
hold on
area(1:N, [shortageCostOverTime; expansionCostOverTime; pumpingCostOverTime]');
legend('Total cost', 'Shortage cost', 'Expansion Cost', 'Pumping Cost')

subplot(2,2,4)
plot(1:N,shortageOverTime)
hold on
plot(1:N,demandOverTime)
plot(1:N,supplyOverTime)
plot(1:N, gwSupplyOverTime)
legend('shortage', 'demand', 'supply', 'gw pumped')

end
%% Save results

if saveOn
    datetime=datestr(now);  
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    mkdir(datetime)
    save(datetime);
end


toc
