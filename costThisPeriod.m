function [ cost, shortageCost, expansionCost, pumpingCost ] = costThisPeriod( a1, a2, costParam, shortage, gw_supply)
% Calculate the cost in the current period given acitons and shortage

% Costs include shortage costs, expansion costs, and pumping costs

shortageCost = shortage * costParam.shortage_cost;
expansionCost = costParam.expansion_cost * a2;
pumpingCost = costParam.pumping_cost * a1 * gw_supply;
cost = shortageCost + expansionCost + pumpingCost;


end

