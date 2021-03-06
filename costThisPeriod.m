function [ cost, shortageCost, expansionCost, pumpingCost, marginalDesalCost ] = costThisPeriod( a1, a2, costParam, shortage, gw_supply, t, s1, exp_supply, gwParam)
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

if gwParam.exaggeratePumpCost
    v= 70;
    cost_per_kwh=.01;
    aquiferDepth = 2;
end

if gwParam.pumpingSubsidy
    cost_per_kwh = .02;
end

pump_cost_perunit=density * 9.81 * (((f*(s1+aquiferDepth)^2*v^2)/(2*9.81*D))+(s1+aquiferDepth)) * conversion_factor ...
    * (100/percent_eff) * cost_per_kwh;

discountFactor = 1/((1+costParam.discount_rate)^t);


% Expansion costs
unitCapex = -0.0014 * water.desal_capacity_expansion.large + 1362.5;
unitOpex = -0.0001 * water.desal_capacity_expansion.large + 186.76

switch a2
    case 0
        expansionCost = 0;
    case 1
        capex = (-0.0014 * water.desal_capacity_expansion.small/365 + 1362.5) * water.desal_capacity_expansion.small/365;
        opex = (-0.0001 * water.desal_capacity_expansion.small/365 + 186.76) * water.desal_capacity_expansion.small/365;
        expansionCost = capex * discountFactor;
    case 2
        capex = (-0.0014 * water.desal_capacity_expansion.large/365 + 1362.5) * water.desal_capacity_expansion.large/365;
        opex = (-0.0001 * water.desal_capacity_expansion.large/365 + 186.76) * water.desal_capacity_expansion.large/365;
        expansionCost = capex * discountFactor;
end
        
% Costs include shortage costs, expansion costs, and pumping costs
shortageCost = shortage * costParam.shortage_cost * discountFactor;
pumpingCost =  pump_cost_perunit * a1 * gw_supply * discountFactor;
marginalDesalCost = exp_supply * costParam.marginal_cost;
cost = shortageCost + expansionCost + pumpingCost + marginalDesalCost;


end

