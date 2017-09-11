function [ cost, shortageCost, expansionCost, pumpingCost ] = costThisPeriod( a1, a2, costParam, shortage, gw_supply, t, s1)
% Calculate the cost in the current period given acitons and shortage

% pumping cost in $/m^3
density=1000; %in kg/m^3
cost_per_kwh=.10;
conversion_factor=2.777778e-7;
percent_eff=80;
f=0.013;
D=0.508;
v=0.24669;
aquiferDepth = 1200;

pump_cost_perunit=density * 9.81 * (((f*(s1+aquiferDepth)^2*v^2)/(2*9.81*D))+(s1+aquiferDepth)) * conversion_factor ...
    * (100/percent_eff) * cost_per_kwh;
 

% Expansion costs
switch a2
    case 0
        expansionCost = 0;
    case 1
        expansionCost = costParam.expansion_cost.capex.small * discountFactor;
    case 2
        expansionCost = costParam.expansion_cost.capex.large * discountFactor;
end
        
% Costs include shortage costs, expansion costs, and pumping costs
discountFactor = 1/((1+costParam.discount_rate)^t);
shortageCost = shortage * costParam.shortage_cost * discountFactor;
pumpingCost =  pump_cost_perunit * a1 * gw_supply * discountFactor;
cost = shortageCost + expansionCost + pumpingCost;


end

