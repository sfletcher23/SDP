%% Simple SDP gw Model

tic
%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 1000;
costParam.expansion_cost = 4000;
costParam.pumping_cost = 1;

% Water paramters
water = struct;
water.demandPerCapita = 2;
water.desal_capacity_initial = 5;
water.desal_capacity_expansion = 20;

% Population parameters
popParam = struct;
popParam.pop_initial = 100;
popParam.min_growth = 0.02;
popParam.max_growth = 0.08;
popParam.max_growth_delta = 0.01;
popParam.discrete_step_pop = 5;
popParam.discrete_step_growth = 0.01;


%% State Definitions for Groundwater and Expansion Decision
% Number states for gw, expansion
s_gw = 1:5;
s_expand = 1:2;
gw_M = length(s_gw); % Groundawter
exp_M = length(s_expand); % Desalination expanded = 2

% Actions: groundwater pumping high or low or none, expand desal
% a1 pumping actions: 0 no pumping, 1 low pumping, 2 high pumping
a_gw_available = [0 1 2];
a_gw_depleted = [0];
% a2 desal actions: 0 no expand, 1 expand
a_expand_available = [0 1];
a_expand_unavailable = [0];


%% Transition Matrix for Groundwater

% Groundwater transition matrix: varies depending on pumping action
% T2: high pumping
T2_gw = zeros(gw_M,gw_M);
T2_gw(1,1) = 1;
for i = 2:gw_M
    T2_gw(i,i) = 0.5;
    T2_gw(i,i-1) = 0.5;
end

% T1: low pumping
T1_gw = zeros(gw_M,gw_M);
T1_gw(1,1) = 1;
for i = 2:gw_M
    T1_gw(i,i) = 0.8;
    T1_gw(i,i-1) = 0.2;
end

% T0: no pumping
T0_gw = zeros(gw_M,gw_M);
T0_gw(1,1) = 1;
for i = 2:gw_M
    T0_gw(i,i) = 1;
    T0_gw(i,i-1) = 0;
end


%% State Definitions and Transition for Pop and Growth

[s_pop, s_growth, T_growth_lookup, nextPop] = gen_pop_growth_states(popParam, N);
pop_M = length(s_pop);
growth_M = length(s_growth);
popParam.growth_initial = s_growth(round(length(s_growth/2)));

%% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = zeros(gw_M, exp_M, pop_M, growth_M, N+1);
X1 = zeros(gw_M, exp_M, pop_M, growth_M, N+1);
X2 = zeros(gw_M, exp_M, pop_M, growth_M, N+1);

% Terminal period
X1(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
X2(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
V(:,:,:,:,N+1) = zeros(gw_M, exp_M, pop_M, growth_M, 1);
    

%% Backwards Recursion

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Loop over all states
    % Loop over groundwater state: 1 is depleted, M1 is full
    for s1 = s_gw       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for s2 = s_expand
            % Loop over population state
            for s3 = s_pop
                % Loop over growth state
                for s4 = s_growth
                
                    bestV = 999999999;
                    bestX = [0 0];  % Groundwater action and expansion action

                    % Update available actions based on whether gw depleted
                    if s1 > 1 
                        a_gw = a_gw_available;
                    else
                        a_gw = a_gw_depleted;
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
                            
                            % Calculate cost this period
                            [shortage, ~] =  shortageThisPeriod(a1, a2, s1, s2, s3, water);
                            cost = costThisPeriod(a1, a2, costParam, shortage);
                            
                            % Caculate state indexes
                            index_s1= find(s1 == s_gw);
                            index_s2 = find(s2 == s_expand);
                            index_s3 = find(s3 == s_pop);
                            index_s4 = find(s4 == s_growth);
                            
                            % Calculate transition matrix
                                
                                % Pick transmat vector for gw based on action, current
                                % gw state
                                if a1 == 0 
                                    T_gw = T0_gw(s1,:);
                                elseif a1 == 1
                                    T_gw = T1_gw(s1,:);
                                else
                                    T_gw = T2_gw(s1,:);
                                end

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
                           expV = sum(T .* nextV);
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
shortageOverTime = zeros(1,N);
supplyOverTime = zeros(1,N);
demandOverTime = zeros(1,N);

% Initial state
s_gw_initial = 5;
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
    [shortageOverTime(t), supplyOverTime(t), demandOverTime(t)] = shortageThisPeriod(action_gw(t), ...
        action_expand(t), state_gw(t), state_expand(t), state_pop(t), water);
    costOverTime(t) = costThisPeriod(action_gw(t), action_expand(t), ...
        costParam, shortageOverTime(t));
    
    % Get transisition mat to next state give current state and actions

        % Get transmat vector to next GW state
        if action_gw(t) == 0
            T_current_gw = T0_gw(state_gw(t),:);
        elseif action_gw(t) == 1
            T_current_gw = T1_gw(state_gw(t),:);
        elseif action_gw(t) == 2
            T_current_gw = T2_gw(state_gw(t),:);
        else
            error('invalid action')
        end
    
        % Get transmat vector for next expansion state (deterministic)
        if action_expand(t) == 1 || state_expand(t) == 2   % desal already expanded or will expand
            T_current_expand = [0 1];
        else
            T__current_expand = [1 0];
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
plot(1:N, state_gw')
hold on
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
legend('Costs')
subplot(3,2,6)
plot(1:N,shortageOverTime)
hold on
plot(1:N,demandOverTime)
plot(1:N,supplyOverTime)
legend('shortage', 'demand', 'supply')

toc