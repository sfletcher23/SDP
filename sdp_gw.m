function [ V, X1, X2, T_gw_all, cumTgw, numRelevantSamples, stateInfeasible, lowestCost, lowestCostAction, s_gw, s_expand, exp_vectors, K_samples, S_samples ] = ...
    sdp_gw( runParam, costParam, popParam, gwParam, water, datetime )
% Run SDP for groundwater model. 


V = [];
X1 = [];
X2 = [];

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
exp_vectors = [];
if runParam.capacityDelay
   s_exp_on = s_expand;
   s_exp_delay1 = s_expand;
   s_exp_delay2 = [s_expand(1) s_expand(2)];
   % Index for feasible expansion state combinations - max total 3v across substates
    s_expand = [1:7 9 10 13 17];
    exp_M = length(s_expand);
    exp_vectors = {s_exp_on' s_exp_delay1' s_exp_delay2'};
end


%% Get K and S samples and use to prune state space

N = runParam.N;
numRelevantSamples = [];
stateInfeasible = [];

if runParam.calculateTgw

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
    y = netscript(x, runParam.adjustOutput);
    minDrawdownHydrograph = y(gwParam.wellIndex,:);
    x = [ones(1,N) * minK; ones(1,N) * minS; time]; 
    y = netscript(x, runParam.adjustOutput);
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
            s1_now = s_gw(index_s1);
            [T_gw_temp, numRel, stateInf, indAbv, indBlw, indRel] = ...
                gw_transrow_nn(gwParam, t, K_samples, S_samples, s1_now, s_gw, runParam.adjustOutput);
            T_gw_all(:,index_s1,t) = T_gw_temp';
            numRelevantSamples(index_s1,t) = numRel;
            stateInfeasible(index_s1,t) = stateInf;
            indexAbove{index_s1, t} = indAbv;
            indexBelow{index_s1, t} = indBlw;
        end
    end    
    
%     save(strcat('T_gw_',datetime), 'T_gw_all', 'K_samples', 'S_samples', 'index_s_gw_time', 'numRelevantSamples', 'stateInfeasible')
    
else
    load(gwParam.TgwLoadName)
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

if runParam.simpleVersion
    
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



end


%% Backwards Recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = NaN(gw_M, exp_M, N+1);
X1 = NaN(gw_M, exp_M, N+1);
X2 = NaN(gw_M, exp_M, N+1);

% Terminal period
X1(:,:,N+1) = zeros(gw_M, exp_M, 1);
X2(:,:,N+1) = zeros(gw_M, exp_M, 1);
V(:,:,N+1) = zeros(gw_M, exp_M, 1);

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
            
            if runParam.capacityDelay
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
            if runParam.capacityDelay
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
                    [ cost, ~,~, ~,~, ~, ~, ~, ~, ~ ] = supplyAndCost( a1, a2, s1, s2, costParam, water, gwParam, t, demandThisPeriod, runParam.capacityDelay, exp_vectors);

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
                    
                    if runParam.capacityDelay
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
            
%             if runParam.saveOn
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


end
%% Solve for optimal policies when all decisions made in 1st stage

% Note: this implementation assumes always pump when you can. Valid for
% normal pumping costs but not exaggerated pumping costs.

lowestCost = [];
lowestCostAction = [];

if runParam.solveNoLearning
    
    % Get water demand
    waterDemand = gwParam.pumpingRate;
    % Get hydrograph for each sample;
    netname = strcat('myNeuralNetworkFunction_', num2str(gwParam.nnNumber));
    netscript = str2func(netname);
    headSample = zeros(length(K_samples), N);
    for i = 1:length(K_samples)
        x = [repmat(K_samples(i),[1,N]); repmat(S_samples(i),[1,N]); [365:365:365*(N)]];
        tempHead = netscript(x, runParam.adjustOutput); 
        headSample(i,:) = tempHead(gwParam.wellIndex,:);
    end
    headSampleRounded = round2x(headSample, s_gw);
    indexPumpingOff = headSample <= gwParam.depthLimit;
    
    meanCostOverTime = zeros(3,N);
    meanShortageOverTime = zeros(3,N);
    meanTotalCost = zeros(3,1);
    meanTotalShortage = zeros(3,1);
    for index_a2 = 1:3
        a2 = a_expand_available(index_a2);
        costOverTime = zeros(gwParam.sampleSize,N);
        expansionCostOverTime = zeros(gwParam.sampleSize,N);
        pumpingCostOverTime = zeros(gwParam.sampleSize,N);
        shortageOverTime = zeros(gwParam.sampleSize,N);
        capacityOverTime = zeros(gwParam.sampleSize,N);
        marginalDesalCostOverTime = zeros(gwParam.sampleSize,N);
        totalCost = zeros(gwParam.sampleSize,1);
        totalShortage = zeros(gwParam.sampleSize,1);
        parfor i = 1:gwParam.sampleSize
            costOverTime_sample = zeros(1,N);
            expansionCostOverTime_sample = zeros(1,N);
            pumpingCostOverTime_sample = zeros(1,N);
            shortageOverTime_sample = zeros(1,N);
            capacityOverTime_sample = zeros(1,N);
            marginalDesalCostOverTime_sample = zeros(1,N);
            for t = 1:N
                s1 = 200 - headSampleRounded(i,t);
                a1 = ~indexPumpingOff(i,t);
                if t ==1
                    s2_now = 0;
                    a2_now = a2;
                else
                    a2_now = 0;
                    switch a2
                        case 0
                            s2_now = 0;
                        case 1
                            s2_now = water.desal_capacity_expansion.small;
                        case 2
                            s2_now = water.desal_capacity_expansion.large;
                    end
                end
                [ cost, shortageCost, expansionCost, pumpingCost, marginalDesalCost, shortage, capacity, minjur_supply, exp_supply, othergw_supply ] ...
                 = supplyAndCost( a1, a2_now, s1, s2_now, costParam, water, gwParam, t, gwParam.pumpingRate, false, []);
                costOverTime_sample(t) = cost;
                expansionCostOverTime_sample(t) = expansionCost;
                pumpingCostOverTime_sample(t) = pumpingCost;
                shortageOverTime_sample(t) = shortage;
                capacityOverTime_sample(t) = capacity;
                marginalDesalCostOverTime_sample(t) = marginalDesalCost;
            end
            costOverTime(i,:) =  costOverTime_sample;
            expansionCostOverTime(i,:) = expansionCostOverTime_sample;
            pumpingCostOverTime(i,:) =  pumpingCostOverTime_sample;
            shortageOverTime(i,:) = shortageOverTime_sample;
            capacityOverTime(i,:) = capacityOverTime_sample;
            marginalDesalCostOverTime(i,:) = marginalDesalCostOverTime_sample;
            totalCost(i)=sum(costOverTime_sample);
            totalShortage(i)=sum(costOverTime_sample);
        end
        meanCostOverTime(index_a2,:) = mean(costOverTime,1);
        meanShortageOverTime(index_a2,:) = mean(shortageOverTime,1);
        meanTotalCost(index_a2) = mean(totalCost);
        meanTotalShortage(index_a2) = mean(totalShortage);
    end
[lowestCost, lowestCostActionIndex] = min(meanTotalCost);
lowestCostAction = a_expand_available(lowestCostActionIndex);
% figure; plot(1:N, meanCostOverTime)
% figure; plot(1:N,  meanShortageOverTime)
    
end




end

