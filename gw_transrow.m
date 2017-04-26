function [ T_gw ] = gw_transrow(s1, s_gw, a1, dd_values, dd_prob, index_s1, index_s3, demand_range, demandThisPeriod)
% Calculates the transition probabilities for the next groundwater state
% given current state, pumping volume (demand), and drawdown values and
% probabilities

T_gw = zeros(1,length(s_gw));

switch a1
    case 1
  
    % Get range of next groundwater values
    next_s = s1 + dd_values;

    % Get index in groundwater state that corresponds to range of possible next
    % values
    [~,next_s_index] = ismember(next_s,s_gw);
    truncate_index = next_s_index == 0;
    next_s_index(truncate_index) = [];

    % Get index in demand_rane that corresponds to current demand
    index_demand = find(abs(demand_range - demandThisPeriod) < 1);
    % check that index_demand and index_s3 are the same
    if index_demand ~= index_s3
        error('demand index and population index are not the same')
    end

    % Get the row from dd_prob for current demand
    dd_prob_row = dd_prob(index_demand,:);

    % Move any probability from values in next_s that are greater than max gw
    % state to the highest possible gw state
    index_max = next_s > max(s_gw);
    max_possible_index = find(index_max,1) - 1;
    if ~isempty(max_possible_index)
        sum_max_prob = sum(dd_prob_row(:,index_max));
        dd_prob_row(index_max) = 0;
        dd_prob_row(max_possible_index) = dd_prob_row(max_possible_index) + sum_max_prob;
        dd_prob_row = dd_prob_row(1:max_possible_index);    
    end

    % Select the row from dd_prob corresponding to the current demand, insert
    % into the relevant states in T_g2
    T_gw(next_s_index) = dd_prob_row;

    case 0
    T_gw(index_s1) = 1;
    otherwise
    error('Invalid a1')
    
    % Test valid prob distribution
    margin = 1E-4;
    err = abs(sum(T_gw) - 1);
    if err > margin
        error('invalid probability distribution for T_gw')
    end
    
end

end



