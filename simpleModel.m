%% Simple model to test MDP Toolbox

% State space only on groundwater level, expansion
% Assume constant population / demand

%% Parameters

% Time period
N = 10;

% Cost paramters
costParam = struct;
costParam.shortage_cost = 25;
costParam.expansion_cost = 100000000; 
costParam.pumping_cost = 10000;
costParam.discount_rate = 0.04;

% Constant population
pop = 4;    % million

numStateVectors = 2;
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
s_gw = 0:10:100;
gw_M = length(s_gw);

% Actions: Pump groundwater this period (full demand) or not
    % 1 is pump, 0 is no pump
gw_A = [0 1]';

% Transition array
T_gw = cell(1,2);
T_gw{2} = diag(ones(1,gw_M)/2);
tempAdd = [zeros(gw_M-1,1), T_gw{2}(2:end, 2:end)];
tempAdd = [tempAdd; zeros(1,gw_M)];
T_gw{2} = T_gw{2} + tempAdd;
T_gw{2}(end,end) = 1;
T_gw{1} = diag(ones(1,gw_M));

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


%% Compute transition matrix

%  T is a 1 x A cell array. 
%  Each cell is a transisiton matrix with a numStateVector dimensional matrix.

T3 = cell(1,sizeA);
for i = 1:sizeA
        T3{i} = zeros(S);
        a = A(i,1);
        b = A(i,2);
        range = gw_M*b+1:gw_M*(b+1);    % select which range is fille depending on exp state
        T3{i}(gw_M+1:gw_M*2,gw_M+1:gw_M*2) = T_gw{a+1};  % 4th quadrant always filled
        T3{i}(1:gw_M, range) = T_gw{a+1}; % Fill relevant 2nd range
        T3{i} = sparse(T3{i});
end


%% Compute reward matrix

% Calculate demand
demandConstant= demand(water, pop, fraction);

% Calculate cost and shortages this period
[shortage, ~, ~, gw_supply] =  shortageThisPeriod(a1, a2, s1, s2, s3, water, demandThisPeriod, s_gw, gwParam);
cost = costThisPeriod(a1, a2, costParam, shortage, gw_supply, t);





