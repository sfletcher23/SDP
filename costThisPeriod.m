function [ cost ] = costThisPeriod( a1, a2, costParam, shortage)
% Calculate the cost in the current period given acitons and shortage

% Costs include shortage costs, expansion costs, and pumping costs
cost = shortage * costParam.shortage_cost + costParam.expansion_cost * a2 ...
    + costParam.pumping_cost * a1;


end

