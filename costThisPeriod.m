function [ cost, shortageCost, expansionCost, pumpingCost ] = costThisPeriod( a1, a2, costParam, shortage, gw_supply, t)
% Calculate the cost in the current period given acitons and shortage

% Costs include shortage costs, expansion costs, and pumping costs
discountFactor = 1/((1+costParam.discount_rate)^t);
shortageCost = shortage * costParam.shortage_cost * discountFactor;
expansionCost = costParam.expansion_cost * a2 * discountFactor;
pumpingCost = costParam.pumping_cost * a1 * gw_supply * discountFactor;
cost = shortageCost + expansionCost + pumpingCost;


end

