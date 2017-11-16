function [ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors )
% Use SDP policy or single time period decision to 


R = runParam.simNum;
N = runParam.N;
[gw_M, exp_M, ~] = size(V);

% Set cost function
if runParam.oldCost
    cost_supply_func = str2func('supplyAndCost_old');
else
    cost_supply_func = str2func('supplyAndCost');
end


% Initialize vector tracking state, actions, water balance, costs over time 
sim = struct;
sim.state_gw = zeros(R,N);
sim.state_expand = zeros(R,N);
sim.action_gw = zeros(R,N);
sim.action_expand = zeros(R,N);
sim.costOverTime = zeros(R,N);
sim.shortageCostOverTime = zeros(R,N);
sim.expansionCostOverTime = zeros(R,N);
sim.pumpingCostOverTime = zeros(R,N);
sim.shortageOverTime = zeros(R,N);
sim.capacityOverTime = zeros(R,N);
sim.minjurSupplyOverTime = zeros(R,N);
sim.othergwSupplyOverTime = zeros(R,N);
sim.demandOverTime = zeros(R,N);
sim.expSupplyOverTime = zeros(R,N);
sim.margDesalCostOverTime = zeros(R,N);
sim.T_gw_time = zeros(gw_M,N,R);
sim.sampleIndexOverTime = zeros(gwParam.sampleSize,N,R);
sim.failureProbOverTime = zeros(R,N);

% Initial state
s_gw_initial = 0;
if runParam.capacityDelay 
    s_expand_initial = 1;
else
    s_expand_initial = 0;
end

state_gw(1) = s_gw_initial;
state_expand(1) = s_expand_initial;

% [K_samples, S_samples] = gen_param_dist(gwParam.infoScenario, gwParam.sampleSize, 1, N);

for i = 1:R
    
    state_gw_now = zeros(1,N);
    state_expand_now = zeros(1,N);
    state_expand_now(1) = s_expand_initial;
    action_gw_now = zeros(1,N);
    action_expand_now = zeros(1,N);
    costOverTime_now = zeros(1,N);
    shortageCostOverTime_now = zeros(1,N);
    expansionCostOverTime_now = zeros(1,N);
    pumpingCostOverTime_now = zeros(1,N);
    shortageOverTime_now = zeros(1,N);
    capacityOverTime_now = zeros(1,N);
    minjurSupplyOverTime_now = zeros(1,N);
    othergwSupplyOverTime_now = zeros(1,N);
    demandOverTime_now = zeros(1,N);
    expSupplyOverTime_now = zeros(1,N);
    margDesalCostOverTime_now = zeros(1,N);
    T_gw_time_now = zeros(gw_M,N);
    sampleIndexOverTime_now = zeros(gwParam.sampleSize,N);
    failureProbOverTime_now = zeros(1,N);
    
    for t = 1:N

        % Caculate state indexes
        index_state_gw = find(state_gw_now(t) == s_gw);
        index_state_expand = find(state_expand_now(t) == s_expand);
        
        % Lookup failure prob if keep pumpin
        failureProbOverTime_now(t) = cumTgw(index_state_gw ,t);
        
        % Lookup optimal policy for current state
        if useNoInfoPolicy
            if t == 1
                action_expand_now(t) = lowestCostAction;
            else 
                action_expand_now(t) = 0;
            end
            if state_expand_now(t) < gwParam.depthLimit && state_expand_now(t) ~= -1
                action_gw_now(t) = 1;
            else
                action_gw_now(t) = 0;
            end
        else
            action_gw_now(t) = X1(index_state_gw, index_state_expand, t);
            action_expand_now(t) = X2(index_state_gw, index_state_expand, t);
        end


        % Calculate demand, shortage, and cost for current t
        demandOverTime_now(t) = gwParam.pumpingRate;
        [ costOverTime_now(t), shortageCostOverTime_now(t), expansionCostOverTime_now(t), pumpingCostOverTime_now(t), margDesalCostOverTime_now(t), ...
            shortageOverTime_now(t), capacityOverTime_now(t),minjurSupplyOverTime_now(t), expSupplyOverTime_now(t), othergwSupplyOverTime_now(t) ] ...
            = cost_supply_func( action_gw_now(t),  action_expand_now(t), state_gw_now(t), state_expand_now(t), ...
            costParam, water, gwParam, t,  demandOverTime_now(t), runParam.capacityDelay, exp_vectors);

        % Get transisition mat to next state give current state and actions

            % Get transmat vector to next GW state 
            if action_gw_now(t) == 0
                T_current_gw = zeros(1, length(s_gw));
                T_current_gw(1) = 1;
                index = ones(gwParam.sampleSize,1)*-99;
            else
                T_current_gw = T_gw_all(:,index_state_gw,t)';
            end

            % Get transmat vector for next expansion state (deterministic)
            T_current_expand = zeros(1,exp_M);
                             
            if runParam.capacityDelay
                subindex_s2 = linIndex2VecIndex(state_expand_now(t), exp_vectors);
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
                if action_expand_now(t) == 1
                    T_exp_delay1_ind = 2;
                elseif action_expand_now(t) == 2
                    T_exp_delay2_ind = 2;
                end
                temp_index = vectorIndex([T_exp_online_ind T_exp_delay1_ind T_exp_delay2_ind], exp_vectors);
                exp_index = find(s_expand == temp_index);
                T_current_expand(exp_index) = 1;
            else
                if action_expand_now(t) == 0
                    T_current_expand(index_state_expand) = 1; % Stay in current state
                elseif action_expand_now(t) == 1
                    T_current_expand(index_state_expand + 1) = 1; % Move up one state
                elseif action_expand_now(t) == 2
                    T_current_expand(index_state_expand + 3) = 1; % Move up three states
                end
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

            state_gw_now(t+1) = s_gw(ind_s1); 
            state_expand_now(t+1) = s_expand(ind_s2);
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
    sim.state_gw(i,:) = state_gw_now;
    sim.state_expand(i,:) = state_expand_now;
    sim.action_gw(i,:) = action_gw_now;
    sim.action_expand(i,:) = action_expand_now;
    sim.costOverTime(i,:) = costOverTime_now;
    sim.shortageCostOverTime(i,:) = shortageCostOverTime_now;
    sim.expansionCostOverTime(i,:) = expansionCostOverTime_now;
    sim.pumpingCostOverTime(i,:) =pumpingCostOverTime_now;
    sim.shortageOverTime(i,:) = shortageOverTime_now;
    sim.capacityOverTime(i,:) = capacityOverTime_now;
    sim.minjurSupplyOverTime(i,:) = minjurSupplyOverTime_now;
    sim.othergwSupplyOverTime(i,:) = othergwSupplyOverTime_now;
    sim.demandOverTime(i,:) = demandOverTime_now;
    sim.expSupplyOverTime(i,:) = expSupplyOverTime_now;
    sim.margDesalCostOverTime(i,:) = margDesalCostOverTime_now;
    sim.T_gw_time(:,:,i) = T_gw_time_now;
    sim.sampleIndexOverTime(:,:,i) = sampleIndexOverTime_now;
    sim.failureProbOverTime(i,:) = failureProbOverTime_now;
end

sim.averageTotalCost = mean(sum(sim.costOverTime,2));
sim.averageTotalShortage = mean(sum(sim.shortageOverTime,2));


end

