%% Simple model to test MDP Toolbox

% State space only on groundwater level, expansion
% Assume constant population / demand

%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 4500;
costParam.expansion_cost = 18000000; 
costParam.pumping_cost = 3025.0;
costParam.discount_rate = 0.01;


numStateVectors = 3;
numActionVectors = 2;


%% Desalination Expansion State vector, transition array, and actions

% Water infrastructure paramters
water = struct;
water.demandPerCapita = 120;    % L/p/d
water.demandFraction = 1/100;
water.desal_capacity_initial = 10E5;
water.desal_capacity_expansion = 5E5;


% State definitions
s_expand = 1:2;
exp_M = length(s_expand); % Desalination expanded = 2

% a2 desal actions: 0 no expand, 1 expand
exp_A = [0 1]';

% Transition array
T_exp = cell(1,2);
T_exp{1} = [1 0 ; 0 1];
T_exp{2} = [0 1 ; 0 1];


%% GW State vector, actions, and transition array

% GW parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.depthLimit = 400;
gwParam.pumpingRate = 7E5;

% States
s_gw = 0:20:gwParam.depthLimit;
gw_M = length(s_gw);

% Actions: Pump groundwater this period (full demand) or not
    % 1 is pump, 0 is no pump
gw_A = [0 1]';

% Transition array
T_gw = cell(1,2);
T_gw{2} = zeros(gw_M);
for i = 1:gw_M
    if i+4 <= gw_M
        T_gw{2}(i,i:i+4) = 1/5;
    elseif i+3 <= gw_M
        T_gw{2}(i,i:i+2) = 1/5;
        T_gw{2}(i,i+3) = 2/5;
    elseif i+2 <= gw_M
        T_gw{2}(i,i:i+1) = 1/5;
        T_gw{2}(i,i+2) = 3/5;
    elseif i+1 <= gw_M
        T_gw{2}(i,i) = 1/5;
        T_gw{2}(i,i+1) = 4/5;
    end
end
T_gw{1} = diag(ones(1,gw_M));

%% Population state vector and transition array

% Population parameters
popParam = struct;
popParam.pop_initial = 4;   % in millions 
popParam.pop_max = 9;
popParam.min_growth = 0.02;
popParam.max_growth = 0.08;
popParam.discrete_step_pop =  0.1;
popParam.discrete_step_growth = 0.01;
popParam.growth_initial = 0.03;

% Get state vector, lenght, and transistion array
[ s_pop, pop_M, T_pop ] = gen_pop_states( popParam);


%% Size of state space
S = gw_M * exp_M * pop_M; 

% Ordered state size
stateSizes = [gw_M exp_M pop_M];


%% Compute action matrix

A_vectors = cell(1,numActionVectors);
A_vectors{1} = gw_A;
A_vectors{2} = exp_A;

% A3 Combine column vectors so that every permutation of the column vectors is
% a row in a matrix with numActionVectors colums
B = A_vectors{1};
[~, col] = size(B);
if col ~= 1
    error('A_vectors must contain only column vectors')
end
for i = 2:length(A_vectors)
    A = A_vectors{i};
    [~, col] = size(A);
    if col ~= 1
        error('A_vectors must contain only column vectors')
    end
    Atemp = repmat(A,1,length(B))';
    [row, col] = size(Atemp);
    Atemp2 = reshape(Atemp,row*col,1);
    [A_length, ~] = size(A);
    Btemp = repmat(B,length(A),1);
    pairs = [Btemp Atemp2];
    B = pairs;
end
A = B;

% Get number of actions in 1D action space
[sizeA, ~] = size(A);
    % Can also use vectorIndex to calculate this
    % vectorIndex can be used to calculate the row index of the combine
    % column vectors based on the indices for each of the individual action
    % vectors


%% Compute transition matrix

%  T is a 1 x A cell array. 
%  Each cell is a transisiton matrix with a numStateVector dimensional matrix.

T = cell(1,sizeA);
T_gw_exp = T;
for i = 1:sizeA
    T_gw_exp{i} = zeros(gw_M * exp_M);
    a = A(i,1);
    b = A(i,2); 
    range = gw_M*b+1:gw_M*(b+1);    % select which range is fille depending on exp state
    T_gw_exp{i}(gw_M+1:gw_M*2,gw_M+1:gw_M*2) = T_gw{a+1};  % 4th quadrant always filled
    T_gw_exp{i}(1:gw_M, range) = T_gw{a+1}; % Fill relevant 2nd range
    %T{i} = sparse(T{i});
end


for i = 1:sizeA
    tempTCell = cell(pop_M);
    for j = 1:pop_M
        for k = 1:pop_M
            tempTCell{j,k} = T_gw_exp{i} * T_pop(j,k);
        end
    end
    T{i} = cell2mat(tempTCell);
end

% Later: generalize this 

%% Compute reward matrix


R = zeros(S, sizeA);

for s1 = s_gw
    index_s1 = find(s1 == s_gw);
    for s2 = s_expand
        index_s2 = find(s2 == s_expand);
        for s3 = s_pop
            index_s3 = find(s3 == s_pop);
            for a1 = gw_A'
                for a2 = exp_A'

                    % If this action is not available at this state, cost = Inf
                    if (s2 == 2 && a2 == 1) || (s1 == s_gw(end) && a1 == 1)
                        cost = 1e30;
                    else % Calculate cost and shortages     
                        demand = get_demand(water, s3, water.demandFraction); % Calculate demand
                        [ cost, ~, ~, ~, ~, ~, ~] = getCost(a1, a2, s1, s2, water, demand, s_gw, gwParam, costParam);  % Can vectorize this fucntion later
                    end

                    % Calculate action index
                    indexAction = vectorIndex([a1+1 a2+1], A_vectors);   % +1 since indexing starts at 1, actions start at 0, maybe change action definitions?

                    % Calculate state vectors
                    indexState = vectorIndex([index_s1 index_s2 index_s3], {s_gw', s_expand', s_pop'});

                    % Put cost in reward matrix
                    R(indexState, indexAction) = -cost;

                end
            end
        end
    end
end


%% Solve Infinite time horizon SDP using value iteration

% Check problem is in standard formulation
error_msg = mdp_check(T , R);

%[policy, iter, cpu_time_value] = mdp_value_iteration (T, R, costParam.discount_rate); 
[V, policy, iter, cpu_time_policy] = mdp_policy_iteration (T, R, costParam.discount_rate); 

%% Visualize results: plot optimal policies

gw_step = s_gw(2) - s_gw(1);
exp_step = s_expand(2) - s_expand(1);
color = {'b', 'g', 'y', 'r'};
fig = figure;

samplePop = [4.3 4.4 4.5 4.6 4.7 4.8];
for i = 1:length(samplePop)
    indexSamplePop(i) = find(samplePop(i) == s_pop);
end
for k = 1:length(samplePop)
    subplot(3,2,k)
    for i = 1:gw_M
        for j = 1:exp_M
            
            if k == length(samplePop)
                patch(1,1,color{1}) % Just to make legend work, will be covered up later
                patch(1,1,color{2})
                patch(1,1,color{3})
                patch(1,1,color{4})
                leg = legend('No pump, no expand', 'Pump, no expand', 'No pump, expand', 'Pump, expand');
%                 leg.Location = 'southeastoutside';
            end
            
            x = [s_gw(i)-(gw_step/2) s_gw(i)-(gw_step/2) s_gw(i)+(gw_step/2) s_gw(i)+(gw_step/2)];
            y = [s_expand(j)-(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)+(exp_step/2) s_expand(j)-(exp_step/2)];
            stateIndex =  vectorIndex([i j indexSamplePop(k)], {s_gw', s_expand', s_pop'});
            policyThisState = policy(stateIndex);
            colorThisState = color{policyThisState};
            patch(x,y,colorThisState)
            hold on  
            
        end
    end
    ax = gca;
    ax.XTick = s_gw;
    ax.YTick = s_expand;
    xlim([s_gw(1)-gw_step/2 s_gw(end)+gw_step/2])
    ylim([s_expand(1)-exp_step/2 s_expand(end)+exp_step/2])
    xlabel('Groundwater state')
    ylabel('Expand state')
    ax.YTickLabel = {'Not expanded', 'Expanded'};
    title(strcat('Population: ', num2str(s_pop(indexSamplePop(k)))))
    ax.XTickLabelRotation = 90;
end


%% Simulate water system
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

N = 20;

% Initialize vector tracking state, actions, water balance, costs over time 
state_gw = zeros(1,N);
state_expand = zeros(1,N);
state_pop = zeros(1,N);
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

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;
state_pop(1) = s_pop_initial;

index_state_gw = find(state_gw(1) == s_gw);
index_state_expand = find(state_expand(1) == s_expand);
index_state_pop = find(state_pop(1) == s_pop);


for t = 1:N
    
    % Lookup optimal policy for current state
    currentState = vectorIndex([index_state_gw index_state_expand index_state_pop], {s_gw', s_expand', s_pop'});
    bestPolicy = policy(currentState);
    action_gw(t) = A(bestPolicy,1);
    action_expand(t) = A(bestPolicy,2);
    
    % Calculate demand, shortage, and cost for current t
    demandOverTime(t) = get_demand( water, state_pop(t), water.demandFraction );
    [costOverTime(t), shortageCostOverTime(t), expansionCostOverTime(t), pumpingCostOverTime(t), ...
        shortageOverTime(t), gwSupplyOverTime(t), supplyOverTime(t)]  = ...
        getCost(action_gw(t), action_expand(t), state_gw(t), state_expand(t), water, demandOverTime(t), s_gw, gwParam, costParam);
        
    % Simulate next state
    if t < N
        p = rand();
        Tcum = cumsum(T{bestPolicy}(currentState,:));
        nextState = find(p < Tcum,1);
        nextStateVectors = linIndex2VecIndex(nextState, {s_gw', s_expand', s_pop'});

        index_state_gw = nextStateVectors(1); 
        index_state_expand = nextStateVectors(2);
        index_state_pop = nextStateVectors(3);
        
        state_gw(t+1) = s_gw(index_state_gw);
        state_expand(t+1) = s_expand(index_state_expand);
        state_pop(t+1) = s_pop(index_state_pop);

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


% Plot system performance
subplot(3,2,5)
plot(1:N,costOverTime);
h = gca;
h.YLim(1) = 0;
hold on
area(1:N, [shortageCostOverTime; expansionCostOverTime; pumpingCostOverTime]');
legend('Total cost', 'Shortage cost', 'Expansion Cost', 'Pumping Cost')

subplot(3,2,4)
plot(1:N,shortageOverTime)
hold on
plot(1:N,demandOverTime)
plot(1:N,supplyOverTime)
plot(1:N, gwSupplyOverTime)
legend('shortage', 'demand', 'supply', 'gw pumped')



