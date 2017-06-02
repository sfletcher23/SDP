%% Simple model to test MDP Toolbox

% State space only on groundwater level, expansion
% Assume constant population / demand

%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 4500;
costParam.expansion_cost = 100000000; 
costParam.pumping_cost = 3001;
costParam.discount_rate = 0.07;


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
gwParam.depthLimit = 100;
gwParam.pumpingRate = 7E5;

% States
s_gw = 0:10:200;
gw_M = length(s_gw);

% Actions: Pump groundwater this period (full demand) or not
    % 1 is pump, 0 is no pump
gw_A = [0 1]';

% Transition array
T_gw = cell(1,2);
T_gw{2} = zeros(gw_M);
for i = 1:gw_M
    if i+3 <= gw_M
        T_gw{2}(i,i:i+3) = 1/4;
    elseif i+2 <= gw_M
        T_gw{2}(i,i:i+1) = 1/4;
        T_gw{2}(i,i+2) = 1/2;
    elseif i+1 <= gw_M
        T_gw{2}(i,i) = 1/4;
        T_gw{2}(i,i+1) = 3/4;
    else
        T_gw{2}(i,i) = 1;
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
                        [ cost, ~, ~, ~ ] = getCost(a1, a2, s1, s2, water, demand, s_gw, gwParam, costParam);  % Can vectorize this fucntion later
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

%% Visualize results

gw_step = s_gw(2) - s_gw(1);
exp_step = s_expand(2) - s_expand(1);
color = {'b', 'g', 'y', 'r'};
fig = figure;

samplePop = [4 4.5 5 5.5 6 6.5];
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
end
