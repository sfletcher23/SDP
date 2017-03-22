%% Simple SDP gw Model

%% Definitions / Initialize

% Time period
N = 20;

% Number states
s_gw = 1:5;
s_expand = 1:2;
M1 = length(s_gw); % Groundawter
M2 = length(s_expand); % Desalination expanded = 2

% Actions: groundwater pumping high or low or none, expand desal
% a1 pumping actions: 0 no pumping, 1 low pumping, 2 high pumping
a_gw_available = [0 1 2];
a_gw_depleted = [0];
% a2 desal actions: 0 no expand, 1 expand
a_expand_available = [0 1];
a_expand_unavailable = [0];

% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = zeros(M1, M2, N+1);
X1 = zeros(M1, M2, N+1);
X2 = zeros(M1, M2, N+1);

% Terminal period
X(:,:,N+1) = zeros(M1, M2, 1);
V(:,:,N+1) = zeros(M1, M2, 1);

% Initial state
s_gw_initial = 5;
s_expand_initial = 1;

%% Parameters

% Cost paramters
costParam = struct;
costParam.shortage_cost = 10;
costParam.expansion_cost = 4000;
costParam.pumping_cost = 5;

% Water paramters
water = struct;
water.demand = 25;
water.desal_capacity_initial = 5;
water.desal_capacity_expansion = 20;



%% Transition Matrix

% T2: high pumping
T2 = zeros(M1,M1);
T2(1,1) = 1;
for i = 2:M1
    T2(i,i) = 0.5;
    T2(i,i-1) = 0.5;
end

% T1: low pumping
T1 = zeros(M1,M1);
T1(1,1) = 1;
for i = 2:M1
    T1(i,i) = 0.8;
    T1(i,i-1) = 0.2;
end

% T0: no pumping
T0 = zeros(M1,M1);
T0(1,1) = 1;
for i = 2:M1
    T0(i,i) = 1;
    T0(i,i-1) = 0;
end


%% Backwards Recursion

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Loop over all states
    % Loop over groundwater state: 1 is depleted, M1 is full
    for s1 = s_gw       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for s2 = s_expand

            bestV = 999999;
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
                    
                    % Calculate shortage this period
                    [shortage, ~] =  shortageThisPeriod(a1, a2, s1, s2, water);
                    
                    % Calculate cost this period
                    cost = costThisPeriod(a1, a2, costParam, shortage);
                    
                    % Pick transmat row for gw based on action, current
                    % gw state
                    if a1 == 0 
                        T_gw = T0(s1,:);
                    elseif a1 == 1
                        T_gw = T1(s1,:);
                    else
                        T_gw = T2(s1,:);
                    end
                    
                    % Choose V based on which desal state
                    if a2 == 1 || s2 == 2   % desal already expanded or will expand
                        nextV_expand = V(:,2,t+1);
                    else
                        nextV_expand = V(:,1,t+1);
                    end
                    
                   % Calculate expected future cost
                   expV = T_gw * nextV_expand;
                                            
                   % Check if best decision
                   checkV = cost + expV;
                   if checkV < bestV
                       bestV = checkV;
                       bestX = [a1 a2];
                   end
                      
                end
            end
            
            % Save best value and action for current state
            V(s1, s2, t) = bestV;
            X1(s1, s2, t) = bestX(1);
            X2(s1, s2, t) = bestX(2);
            
        end
    end
end

% Why is pumping the same no matter whether expand or not expand?
% Why does desal expand no matter what groundwater state?

%% Simulate performance
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

% Initialize vector tracking state, actions, water balance, costs over time 
state_gw = zeros(1,N);
state_expand = zeros(1,N);
action_gw = zeros(1,N);
action_expand = zeros(1,N);
costOverTime = zeros(1,N);
shortageOverTime = zeros(1,N);
supplyOverTime = zeros(1,N);

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;


for t = 1:N
    
    % Lookup optimal policy for current state
    action_gw(t) = X1(state_gw(t), state_expand(t), t);
    action_expand(t) = X2(state_gw(t), state_expand(t), t);
    
    % Get transition p to next GW state
    if action_gw(t) == 0
        T_current_gw = T0(state_gw(t),:);
    elseif action_gw(t) == 1
        T_current_gw = T1(state_gw(t),:);
    elseif action_gw(t) == 2
        T_current_gw = T2(state_gw(t),:);
    else
        error('invalid action')
    end
    
    % Calculate shortage and cost
    [shortageOverTime(t), supplyOverTime(t)] = shortageThisPeriod(action_gw(t), action_expand(t), ...
        state_gw(t), state_expand(t),water);
    costOverTime(t) = costThisPeriod(action_gw(t), action_expand(t), ...
        costParam, shortageOverTime(t));
    
    if t < N
    % Simulate next GW state
    p = rand();
    p_in_range = false;
    test_state = 1;
    cumProb = 0;
    while ~p_in_range
        cumProbNext = cumProb + T_current_gw(test_state);
        if p > cumProb && p <= cumProbNext
            next_state = test_state;
            break
        end
        test_state = test_state + 1;
        cumProb = cumProbNext;
        if test_state > length(s_gw)
            error('exceeded available states')
        end
    end
    state_gw(t+1) = test_state;
    
    % Next expand state
    if state_expand(t) == 1 && action_expand(t) == 0
        state_expand(t+1) = 1;
    else
        state_expand(t+1) = 2;
    end   
    end
    
end

%% Plot simulation results

figure;
subplot(2,2,1)
plot(1:N, state_gw')
hold on
scatter(1:N, action_gw)
xlabel('time')
legend('GW state', 'pumping level')
subplot(2,2,2)
plot(1:N, state_expand')
hold on
scatter(1:N, action_expand')
xlabel('time')
legend('Expansion state', 'Expansion decision')
subplot(2,2,3)
plot(1:N,costOverTime);
h = gca;
h.YLim(1) = 0;
legend('Costs')
subplot(2,2,4)
plot(1:N,shortageOverTime)
hold on
plot(1:N,water.demand*ones(1,N))
plot(1:N,supplyOverTime)
legend('shortage', 'demand', 'supply')

