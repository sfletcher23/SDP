%% Screening model

%% Screening model inputs and test actions

% Low, medium, or high demand;
demandLevel = 'low';
headDepletion = 'high';

% Test actions
action_gw = ones(1,N);
action_expand = zeros(1,N);
action_expand(1) = 0;
    % Check expansion action valid
    sumExp = sum(action_expand);
    if ~ (sumExp == 0 || sumExp == 1)
        error('invalid expanion option')
    end
    
    
    
    
%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 2;
costParam.expansion_cost = 1e8; 
costParam.pumping_cost = 1;
costParam.discount_rate = 0.04;

% Water infrastructure paramters
water = struct;
water.demandPerCapita = 120;    % L/p/d
water.demandFraction = 1/150;
water.desal_capacity_initial = 10E5;
water.desal_capacity_expansion = 5E5;

% Population parameters
popParam = struct;
popParam.pop_initial = 4;   % in millions 
popParam.min_growth = 0.02;
popParam.max_growth = 0.08;
popParam.max_growth_delta = 0.01;
popParam.discrete_step_pop =  0.07;
popParam.discrete_step_growth = 0.005;
popParam.growth_initial = 0.03;


% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 10000;
gwParam.depthLimit = 5000;
gwParam.pumpingRate = 7E5;

%% State Definitions

% Generate population growth states based on discreitzation above, check
% that the error is less than 1/2 discretization given.
[s_pop, s_growth, T_growth_lookup, nextPop, max_end_pop] = gen_pop_growth_states(popParam, N);
pop_M = length(s_pop);
growth_M = length(s_growth);

% For now, assume some percentage of demand per capita comes from single well
fraction = water.demandFraction;
demand_range = demand(water, s_pop, fraction); % in m^3/y

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
    % 1 is pump, 0 is no pump
a_gw_available = [0 1];
a_gw_depleted = [0];


%% Setup state variables based on screening assumptions

% Initialize vector tracking state, water balance, costs over time 
state_gw = zeros(1,N);
state_expand = zeros(1,N);
state_pop = zeros(1,N);
state_growth = zeros(1,N);
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

% Population and growth states
if strcmp(demandLevel, 'low')
    state_growth(:) = min(s_pop);
elseif strcmp(demandLevel, 'medium')
    t_change = floor(N/2);
    state_growth(1:t_change) = min(s_pop);
    state_growth(t_change+1:end) = max(s_pop);
elseif strcmp(demandLevel, 'high')
    state_growth(:) = max(s_pop);
end
for t = 2:N
    state_pop(t) = round2x( state_pop(t-1) * (1+state_growth(t-1)), s_pop);
end

% Expansion state
state_expand = ones(1,N);
expInd = find(action_expand == 1);
if ~isempty(expInd)
    state_expand(expInd+1:end) = 2;
end

% Groundwater state
for t = 2:N
    ind = find(s_gw == state_gw(t-1));
    if action_gw(t-1) == 1
        if strcmp(headDepletion, 'high')
            state_gw(t) = s_gw(ind+2);
        elseif strcmp(headDepleteion, 'medium')
            state_gw(t) = s_gw(ind+1);
        elseif strcmp(headDepletion, 'low')
            state_gw(t) = s_gw(ind);
        else
            error('Invalid headDepletion value')
        end
    elseif action_gw(t-1) == 0 
        state_gw(t) = state_gw(t-1);
    else
        error('Invalid pumping action')
    end
end
        
           
%% Simulate performance

for t = 1:N
    
    % Caculate state indexes
    index_state_gw = find(state_gw(t) == s_gw);
    index_state_expand = find(state_expand(t) == s_expand);
    index_state_pop = find(state_pop(t) == s_pop);
    index_state_growth = find(state_growth(t) == s_growth);
    
    % Calculate demand, shortage, and cost for current t
    demandOverTime(t) = demand( water, state_pop(t), fraction);
    [shortageOverTime(t), supplyOverTime(t), ~, gwSupplyOverTime(t)] = shortageThisPeriod(action_gw(t), ...
        action_expand(t), state_gw(t), state_expand(t), state_pop(t), water, demandOverTime(t), s_gw, gwParam);
    [costOverTime(t), shortageCostOverTime(t), expansionCostOverTime(t), pumpingCostOverTime(t)]  = ...
        costThisPeriod(action_gw(t), action_expand(t), costParam, shortageOverTime(t),gwSupplyOverTime(t),t);

end

totalCost = sum(costOverTime(t))

%% Plot simulation results

% Plot state evolution w/ actions
figure;

subplot(3,2,1)
% yyaxis left
plot(1:N, state_gw')
% yyaxis right
hold on
scatter(1:N, action_gw*10)
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

