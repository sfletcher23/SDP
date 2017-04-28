%% Simple SDP gw Model

tic

%test comment
%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 10000;
costParam.expansion_cost = 1000000000000; 
costParam.pumping_cost = 10000;

% Water paramters
water = struct;
water.demandPerCapita = 120;    % L/p/d
water.desal_capacity_initial = 7E5;
water.desal_capacity_expansion = 5E5;

% Population parameters
popParam = struct;
popParam.pop_initial = 4;   % in millions 
popParam.min_growth = 0.02;
popParam.max_growth = 0.08;
popParam.max_growth_delta = 0.01;
popParam.discrete_step_pop = 0.1;
popParam.discrete_step_growth = 0.005;

% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 10000;
gwParam.depthLimit = 5000;
gwParam.pumpingRate = 5E5;



%% State Definitions and Transition for Pop and Growth

[s_pop, s_growth, T_growth_lookup, nextPop, max_end_pop] = gen_pop_growth_states(popParam, N);
pop_M = length(s_pop);
growth_M = length(s_growth);
popParam.growth_initial = 0.03;


%% Calculate Demand

% For now, assume some percentage of demand per capita comes from single
% well
fraction = 1/100;
demand_range = demand(water, s_pop, fraction); % in m^3/y

%% State and Action Definitions for Groundwater 


% Import groundwater well data and aquifer properties
% T in units m^2/sec, drawdown and depth in unit meters
groundwaterWells = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A6:I19', 'ReadVariableNames', true, 'ReadRowNames', true);
aquifer = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A25:J27', 'ReadVariableNames', true, 'ReadRowNames', true);

% Generate state space for groundwater head and demand range
[s_gw, gw_M, drawdownMaxAnnual, gwParam.stepSize] = gen_water_growth_states(gwParam, N, demand_range, water, ...
    groundwaterWells, aquifer);

% Actions: Pump groundwater this period (full demand) or not
a_gw_available = [0 1];
a_gw_depleted = [0];


%% Calculate kernel functions
% Calculate kernel functions measuring the response to a unit impulse of
% pumping from each well over time. 
% Kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix

[kernel, T_S_pairs] = gen_kernel(groundwaterWells, aquifer);


%% Calculate pt(k), generate samples of k for each t

index_T_S_samples = zeros(gwParam.sampleSize,N);
for t = 1:N
    % Define pt(k)
        [T_S_pair_cdf] = gen_param_dist(T_S_pairs, t);
    % Generate samples from pt(K)
        p = rand(1,gwParam.sampleSize);
        index_T_S_samples(:,t) = arrayfun(@(x) find(x < T_S_pair_cdf,1), p);
        clear T_S_pair_cdf p
end
%% Desalination Expansions: State Definitions and Actions

% State definitions
s_expand = 1:2;
exp_M = length(s_expand); % Desalination expanded = 2

% a2 desal actions: 0 no expand, 1 expand
a_expand_available = [0 1];
a_expand_unavailable = [0];


%% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = NaN(gw_M, exp_M, pop_M, growth_M, N+1);
X1 = NaN(gw_M, exp_M, pop_M, growth_M, N+1);
X2 = NaN(gw_M, exp_M, pop_M, growth_M, N+1);

% Terminal period
X1(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
X2(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
V(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
    

%% Backwards Recursion

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Calculate range of possible population states this time step
    s_pop_thisPeriod = pop_states_this_period(s_pop, t, nextPop);
        % Check that subset of total pop states
        isSubset = ismembertol(s_pop_thisPeriod,s_pop, 1E-4);
        test = sum(~isSubset);
        if test > 0
            error('Pop growth state set this period invalid')
        end
    
    % Loop over all states
    % Loop over groundwater state: 1 is depleted, M1 is full
    for s1 = s_gw       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for s2 = s_expand
            % Loop over population state
            for s3 = s_pop_thisPeriod
                % Loop over growth state
                for s4 = s_growth
                
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

                    % Loop over all actions
                    % Loop over groundwater pumping action
                    for a1 = a_gw
                        % Loop over expansion action: 1 is do not expand, 2 is expand
                        for a2 = a_expand
                            
                            % Calculate demand
                            demandThisPeriod = demand(water, s3, fraction);
                            
                            % Calculate cost and shortages this period
                            [shortage, ~, ~, gw_supply] =  shortageThisPeriod(a1, a2, s1, s2, s3, water, demandThisPeriod, s_gw, gwParam);
                            cost = costThisPeriod(a1, a2, costParam, shortage, gw_supply);
                            
                            % Caculate state indexes
                            index_s1= find(s1 == s_gw);
                            index_s2 = find(s2 == s_expand);
                            index_s3 = find(abs(s3 - s_pop) < 1E-3);
                            index_s4 = find(s4 == s_growth);
                            
                            % Calculate transition matrix
                                
                                % Get transmat vector for gw based on action, current
                                % gw state
                                T_gw = gw_transrow_kernel(gw_supply, kernel, index_T_S_samples(:,t), t, s1, s_gw );

                                % Get transmat vector for next expansion state
                                % (deterministic)
                                if a2 == 1 || s2 == 2   % desal already expanded or will expand
                                    T_expand = [0 1];
                                else
                                    T_expand = [1 0];
                                end
                                
                                % Get transmat vector for next growth state
                                T_growth = T_growth_lookup(index_s4,:);
                                
                                % Get transmat vector for next population state
                                T_pop = zeros(1, pop_M);
                                nextPopCurrent = nextPop(index_s3, index_s4);   % Value of next pop given current state
                                index_T = find(s_pop == nextPopCurrent);
                                T_pop(index_T) = 1;
                                
                                % Calculate full transition matrix
                                % T gives probability of next state given
                                % current state and actions
                                
                                TRows = cell(4,1);
                                TRows{1} = T_gw;
                                TRows{2} = T_expand;
                                TRows{3} = T_pop;
                                TRows{4} = T_growth;
                                [ T ] = transrow2mat( TRows );
                             
                             % Calculate expected future cost
                           nextV = V(:,:,:,:,t+1);
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

                        end
                    end
                    
                    % Check that bestV is not Inf
                    if bestV == Inf
                        error('BestV is Inf, did not pick an action')
                    end
                    
                    % Save best value and action for current state
                    V(index_s1, index_s2, index_s3, index_s4, t) = bestV;
                    X1(index_s1, index_s2, index_s3, index_s4, t) = bestX(1);
                    X2(index_s1, index_s2, index_s3, index_s4, t) = bestX(2);

                end
            end
        end
    end
end



%% Simulate performance
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

% Initialize vector tracking state, actions, water balance, costs over time 
state_gw = zeros(1,N);
state_expand = zeros(1,N);
state_pop = zeros(1,N);
state_growth = zeros(1,N);
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
s_pop_initial = popParam.pop_initial;
s_growth_initial = popParam.growth_initial;

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;
state_pop(1) = s_pop_initial;
state_growth(1) = s_growth_initial;

for t = 1:N
    
    % Caculate state indexes
    index_state_gw = find(state_gw(t) == s_gw);
    index_state_expand = find(state_expand(t) == s_expand);
    index_state_pop = find(state_pop(t) == s_pop);
    index_state_growth = find(state_growth(t) == s_growth);
    
    % Lookup optimal policy for current state
    action_gw(t) = X1(index_state_gw, index_state_expand, index_state_pop, index_state_growth, t);
    action_expand(t) = X2(index_state_gw, index_state_expand, index_state_pop, index_state_growth, t);
    
    % Calculate demand, shortage, and cost for current t
    demandOverTime(t) = demand( water, state_pop(t), fraction);
    [shortageOverTime(t), supplyOverTime(t), ~, gwSupplyOverTime(t)] = shortageThisPeriod(action_gw(t), ...
        action_expand(t), state_gw(t), state_expand(t), state_pop(t), water, demandOverTime(t), s_gw, gwParam);
    [costOverTime(t), shortageCostOverTime(t), expansionCostOverTime(t), pumpingCostOverTime(t)]  = ...
        costThisPeriod(action_gw(t), action_expand(t), costParam, shortageOverTime(t),gwSupplyOverTime(t));
    
    % Get transisition mat to next state give current state and actions

        % Get transmat vector to next GW state 
        T_current_gw = gw_transrow_kernel(gwSupplyOverTime(t), kernel, index_T_S_samples(:,t), t, state_gw(t), s_gw );
 
        % Get transmat vector for next expansion state (deterministic)
        if action_expand(t) == 1 || state_expand(t) == 2   % desal already expanded or will expand
            T_current_expand = [0 1];
        else
            T_current_expand = [1 0];
        end
        
        % Get transmat vector for next growth state
        T_current_growth = T_growth_lookup(index_state_growth,:);

        % Get transmat vector for next population state
        T_current_pop = zeros(1, pop_M);
        nextPopCurrent = nextPop(index_state_pop, index_state_growth);   % Value of next pop given current state
        index_T = find(s_pop == nextPopCurrent);
        T_current_pop(index_T) = 1;
        
        % Get Transition Matrix from rows
        TRows_current = cell(4,1);
        TRows_current{1} = T_current_gw;
        TRows_current{2} = T_current_expand;
        TRows_current{3} = T_current_pop;
        TRows_current{4} = T_current_growth;
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
        state_pop(t+1) = s_pop(ind_s3);
        state_growth(t+1) = s_growth(ind_s4);
            % Test next state
            test_gw = T_current_gw(ind_s1) >= 0;
            test_expand = T_current_expand(ind_s2) >= 0;
            test_pop = T_current_pop(ind_s3) >= 0;
            test_growth = T_current_growth(ind_s4) >= 0; 
            if ~test_gw
                error('Invalid gw state tranisition')
            end
            if ~test_expand
                error('Invalid expand state tranisition')
            end
            if ~test_pop
                error('Invalid pop state tranisition')
            end
            if ~test_growth
                error('Invalid growth state tranisition')
            end
    end
    
end

%% Plot simulation results

% Plot state evolution w/ actions
figure;
subplot(3,2,1)
yyaxis left
plot(1:N, state_gw')
hold on
yyaxis right
scatter(1:N, action_gw)
xlabel('time')
legend('GW state', 'pumping level')

subplot(3,2,2)
plot(1:N, state_expand')
hold on
scatter(1:N, action_expand')
xlabel('time')
legend('Expansion state', 'Expansion decision')

subplot(3,2,3)
plot(1:N, state_pop')
legend('Population state')

subplot(3,2,4)
plot(1:N, state_growth')
legend('Growth state')


% Plot system performance
subplot(3,2,5)
plot(1:N,costOverTime);
h = gca;
h.YLim(1) = 0;
hold on
area(1:N, [shortageCostOverTime; expansionCostOverTime; pumpingCostOverTime]');
legend('Total cost', 'Shortage cost', 'Expansion Cost', 'Pumping Cost')

subplot(3,2,6)
plot(1:N,shortageOverTime)
hold on
plot(1:N,demandOverTime)
plot(1:N,supplyOverTime)
plot(1:N, gwSupplyOverTime)
legend('shortage', 'demand', 'supply', 'gw pumped')

% Plot parameter changes over time
figure;
for i = 1:N
    Tsamples = T_S_pairs(index_T_S_samples(:,i),1);
    Ssamples = T_S_pairs(index_T_S_samples(:,i),2);
    hold on
    subplot(1,2,1)
    cdfplot(Tsamples)
    hold on
    subplot(1,2,2)
    cdfplot(Ssamples)
end
leg = arrayfun(@num2str, 1:N, 'UniformOutput', false);
subplot(1,2,1)
title('T distribution over time')
legend(leg)
subplot(1,2,2)
title('S distribution over time')
legend(leg)


%%
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore
mkdir(datetime)
save(datetime);


toc
