function [ cost, shortageCost, expansionCost, pumpingCost, shortage, gw_supply, supply ]...
    = getCost(a1, a2, s1, s2, water, demand, s_gw, gwParam, costParam)
% Calculates shortage in current period gives states, actions, and demand

% Groundwater supplied this period
if a1 == 0 || s1 == max(s_gw)
    gw_supply = 0;
else
    gw_supply = gwParam.pumpingRate;
end
	
% Total water supplied is sum of groundwater and desalinated water
supply = gw_supply + water.desal_capacity_initial + ...
    water.desal_capacity_expansion * (s2-1); % desal if available
 
% Calculate shortage
shortage = max(0, demand - supply);

% Costs include shortage costs, expansion costs, and pumping costs
shortageCost = shortage * costParam.shortage_cost;
expansionCost = costParam.expansion_cost * a2;
pumpingCost = costParam.pumping_cost * a1 * gw_supply;
cost = shortageCost + expansionCost + pumpingCost;




end

% update this for gw supply to reflect withdawal demand insted of a1
 
