%% Simple model to test MDP Toolbox

% State space only on groundwater level, expansion
% Assume constant population / demand

tic
%% Parameters

% Control parameters
readPopStates = true; % Read pop states from spreadsheet?
simulationModelOn = false;
plotsOn = true;
calculateKernels = true;
useSimpleGwT = false;

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 20000;
costParam.expansion_cost = 100000000; 
costParam.pumping_cost = 1200.0;
costParam.discount_rate = 0.07;


numStateVectors = 3;


%% Desalination Expansion State vector, transition array, and actions

% Water infrastructure paramters
water = struct;
water.demandPerCapita = 120;    % L/p/d
water.demandFraction = 1/200;
water.desal_capacity_initial = 10E5;
water.desal_capacity_expansion = 5E5;


% State definitions
s_expand = (1:2)';
exp_M = length(s_expand); % Desalination expanded = 2

% a2 desal actions: 0 no expand, 1 expand
exp_A = [0 1]';

% Transition array
T_exp = cell(1,2);
T_exp{1} = [1 0 ; 0 1];
T_exp{2} = [0 1 ; 0 1];

%% Population state vector and transition array

if readPopStates
     [s_pop, pop_M, T_pop] = gen_pop_states_fromxls(); 
    
else
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

end

%% GW State vector, actions, and transition array

% Assume multiple identical pumping wells

% GW parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.depthLimit = 400;
gwParam.pumpingRate = 7E6;
gwParam.stepSize = 5;
gwParam.sampleSize = 1000; 
gwParam.numPumpWells = 2;

% Actions: Pump groundwater this period (full demand) or not
    % 1 is pump, 0 is no pump
gw_A = [0 1]';


% If not using kernel funcitons, calculate a simple representative
% transition matrix for groundwater
if useSimpleGwT
    % States
    s_gw = (0:gwParam.stepSize:gwParam.depthLimit)';
    gw_M = length(s_gw);

    % Transition array
    T_gw = cell(1,2);
    T_gw{2} = sparse(gw_M,gw_M);
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
    T_gw{1} = spdiags(ones(gw_M,1),0,gw_M,gw_M);
    gw_M_well = gw_M;

% Otherwise, calculate groundwater T from kernels    
else
    % Calculate new kernels if not using existing kernels (for now, assume
    % same kernsl for each well, no interaction effects)
    if calculateKernels
        
        % Calculate kernels and samples
        [s_gw, gw_M, kernel, index_T_S_samples] = waterStatesKernels(gwParam, s_pop, water);
        
        % Calculate T
        T_gw = cell(1,2);
        T_gw{1} = spdiags(ones(gw_M,1),0,gw_M,gw_M);
        temp = gw_transmat_kernel(gwParam.pumpingRate, kernel, index_T_S_samples, s_gw, gw_M );
        T_gw{2} = sparse(temp);
        
        
    end
end

% Update s_gw, gw_M, adn T_gw to include multiple wells

if gwParam.numPumpWells > 1
    
    % Calculate new s_gw and gw_M
    s_gw_new = s_gw;
    gw_M_new = gw_M;
    for i = 2:gwParam.numPumpWells
        temp = repmat(s_gw_new, [gw_M 1]);
        temp2 = reshape(repmat(s_gw', [gw_M_new, 1]), [gw_M * gw_M_new, 1]);
        s_gw_new = [temp temp2];
        [gw_M_new, ~] = size(s_gw_new);
    end
    s_gw = s_gw_new;
    gw_M_well = gw_M;
    gw_M = gw_M_new;
    clear gw_M_new s_gw_new temp temp2 
    
    % Calculate new T_gw for 2 wells
    
    [T_gw] = gw_transmat_multiwell(T_gw, gw_M, gw_M_well);
    
end


%% Size of state space
S = gw_M * exp_M * pop_M; 

% Ordered state size
stateSizes = [gw_M exp_M pop_M];


%% Compute action matrix

A_vectors = cell(1,1);
for i = 1:gwParam.numPumpWells
    A_vectors{i} = gw_A;
end
A_vectors{end+1} = exp_A;

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
    T_gw_exp{i} = sparse(gw_M * exp_M, gw_M * exp_M);
    a = A(i,1);
    b = A(i,2); 
    range = gw_M*b+1:gw_M*(b+1);    % select which range is fille depending on exp state
    T_gw_exp{i}(gw_M+1:gw_M*2,gw_M+1:gw_M*2) = T_gw{a+1};  % 4th quadrant always filled
    T_gw_exp{i}(1:gw_M, range) = T_gw{a+1}; % Fill relevant 2nd range
end

% clear T_gw T_exp

for i = 1:sizeA
    tempTCell = cell(pop_M);
    for j = 1:pop_M
        for k = 1:pop_M
            tempTCell{j,k} = T_gw_exp{i} * T_pop(j,k);
        end
    end
    T{i} = cell2mat(tempTCell);
end

% Later: test this 

% clear T_gw_exp T_pop

%% Compute reward matrix

R = zeros(S, sizeA);

% Calculate demand
% Population based on S index
indexPop = ceil((1:S) ./ (gw_M * exp_M));
pop = s_pop(indexPop,1);
dmd = water.demandPerCapita * 365/1000  ...    % m^3/p/year
    * pop * 1E6 ...   % p
    * water.demandFraction;

% Calculate supply
% Initial desal supply
supply = ones(size(R)) * water.desal_capacity_initial;
% Add expansion supply
ind1 = mod(1:S,gw_M*2) > gw_M;
ind2 = mod(1:S,gw_M*2) == 0;
indexExpanded = ind1 | ind2;
supply(indexExpanded,:) = supply(indexExpanded,:) + water.desal_capacity_expansion; 
% Add pumping supply
numPumpWells = sum(A(:,1:gwParam.numPumpWells),2)';
gw_supply = repmat(numPumpWells,[S,1]) * gwParam.pumpingRate;
supply = supply + gw_supply;

% Calculate shortage
shortage = max(0, dmd - supply);

% Costs include shortage costs, expansion costs, and pumping costs
shortageCost = shortage * costParam.shortage_cost;
expansionCost = zeros(size(R));
indexExp = [A(:,end) == 1]';
expansionCost(:,indexExp) = costParam.expansion_cost;
pumpingCost = costParam.pumping_cost * gw_supply;
R = -(shortageCost + expansionCost + pumpingCost);

% Infeasible actions
infeasibleCost = 1e30;
% Pumping with max drawdown
for i = 1:gwParam.numPumpWells
    indexPerGWM = s_gw(:,i) == max(s_gw(:,i));
    indexMaxDrawdown = repmat(indexPerGWM, [S/length(indexPerGWM), 1]);
    indexPump = A(:,i) == 1;
    R(indexMaxDrawdown, indexPump) = -infeasibleCost;    
end
% Expanding when already expanded
R(indexExpanded,indexExp) = -infeasibleCost;

%% Solve Infinite time horizon SDP using value iteration

% Check problem is in standard formulation
error_msg = mdp_check(T , R);

% Turn on verbose 
%mdp_verbose()
mdp_silent()

%[policy, iter, cpu_time_value] = mdp_value_iteration (T, R, costParam.discount_rate); 
[V, policy, iter, cpu_time_policy] = mdp_policy_iteration (T, R, costParam.discount_rate);
%[V_mod, policy_mod, iter_mod, cpu_time_policy_modified] = mdp_policy_iteration_modified(T, R, costParam.discount_rate, 1e-5); this does not seem to work well
%[Q, V, policy, mean_discrepancy] = mdp_Q_learning(T, R, costParam.discount_rate); took a really long time for 2 wells so canceled it
%[policy, iter, cpu_time] = mdp_value_iterationGS(T, R, costParam.discount_rate);


%% Visualize results: plot optimal policies

if plotsOn
        gw_step = s_gw(2,1) - s_gw(1,1);
        exp_step = s_expand(2) - s_expand(1);
        color = {'b', 'g', 'y', 'r'};
        fig = figure;

        samplePop = [6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11];
        for i = 1:length(samplePop)
            indexSamplePop(i) = find(samplePop(i) == s_pop(:,1),1);
        end
        for k = 1:length(samplePop)
            subplot(4,3,k)
            for i = 1:gw_M_well
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
                    stateIndex =  vectorIndex([i j indexSamplePop(k)], {s_gw, s_expand, s_pop});
                    policyThisState = policy(stateIndex);
                    colorThisState = color{policyThisState};
                    patch(x,y,colorThisState)
                    hold on  

                end
            end
            ax = gca;
            ax.XTick = s_gw(1:gw_M_well,1);
            ax.YTick = s_expand;
            xlim([s_gw(1)-gw_step/2 s_gw(end)+gw_step/2])
            ylim([s_expand(1)-exp_step/2 s_expand(end)+exp_step/2])
            xlabel('Groundwater state')
            ylabel('Expand state')
            ax.YTickLabel = {'Not expanded', 'Expanded'};
            title(strcat('Population: ', num2str(s_pop(indexSamplePop(k)))))
            ax.XTickLabelRotation = 90;
        end
end
%% Simulate water system
% SDP above finds optimal policy for each state and time period. Now, use
% intial state, and transition matrix to simulate performance of the
% system

if simulationModelOn && plotsOn

    N = 20;

    % Initialize vector tracking state, actions, water balance, costs over time 
    state_gw = zeros(1,N);
    state_expand = zeros(1,N);
    state_pop = zeros(2,N);
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
    s_pop_initial = s_pop(1,:);

    state_gw(1) = s_gw_initial;
    state_expand(1) = s_expand_initial;
    state_pop(:,1) = s_pop_initial';

    index_state_gw = find(state_gw(1) == s_gw);
    index_state_expand = find(state_expand(1) == s_expand);
    index_state_pop = find((state_pop(1,1) == s_pop(:,1)) .* (state_pop(2,1) == s_pop(:,2)));


    for t = 1:N

        % Lookup optimal policy for current state
        currentState = vectorIndex([index_state_gw index_state_expand index_state_pop], {s_gw, s_expand, s_pop});
        bestPolicy = policy(currentState);
        action_gw(t) = A(bestPolicy,1);
        action_expand(t) = A(bestPolicy,2);

        % Calculate demand, shortage, and cost for current t
        demandOverTime(t) = get_demand( water, state_pop(1,t), water.demandFraction );
        [costOverTime(t), shortageCostOverTime(t), expansionCostOverTime(t), pumpingCostOverTime(t), ...
            shortageOverTime(t), gwSupplyOverTime(t), supplyOverTime(t)]  = ...
            getCost(action_gw(t), action_expand(t), state_gw(t), state_expand(t), water, demandOverTime(t), s_gw, gwParam, costParam);

        % Simulate next state
        if t < N
            p = rand();
            Tcum = cumsum(T{bestPolicy}(currentState,:));
            nextState = find(p < Tcum,1);
            nextStateVectors = linIndex2VecIndex(nextState, {s_gw, s_expand, s_pop});

            index_state_gw = nextStateVectors(1); 
            index_state_expand = nextStateVectors(2);
            index_state_pop = nextStateVectors(3);

            state_gw(t+1) = s_gw(index_state_gw);
            state_expand(t+1) = s_expand(index_state_expand);
            state_pop(:,t+1) = s_pop(index_state_pop,:)';

        end

    end

end
%% Plot simulation results

if simulationModelOn
    
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
    yyaxis left
    plot(1:N, state_pop(1,:)')
    hold on
    yyaxis right
    plot(1:N, state_pop(2,:)')
    legend('Population state', 'Growth rate state')


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

end

S
runTime = toc
cpu_time_policy
nonzeroel = nnz(T{4})
whos T
