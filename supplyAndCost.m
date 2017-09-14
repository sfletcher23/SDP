function [ cost, shortageCost, expansionCost, pumpingCost, marginalDesalCost, shortage, capacity, minjur_supply, exp_supply, othergw_supply ] ...
    = supplyAndCost( a1, a2, s1, s2, costParam, water, gwParam, t, demand)


% Calculate the capacity, supply, and cost in the current period given
% acitons and state

% First: Calculate capacity available

    % Groundwater capacity this period
    othergw_supply = gwParam.otherPumpingRate;
    if a1 == 0 || s1 == 200 || s1 == -1
        minjur_capacity = 0;
    else
        minjur_capacity = gwParam.pumpingRate;
    end

    % Exp desal capacity
    exp_capacity = s2; % assume use 100% of desal capacity

    % Total water capacity is sum of groundwater and desalinated water
    % capacity = minjur_supply + water.desal_capacity_initial + exp_capacity + othergw_supply;
    capacity = minjur_capacity + exp_capacity ; 

% Second: Calculate expansion costs and marginal costs

    % Discount factor
    discountFactor = 1/((1+costParam.discount_rate)^t);

    % Pumping cost in $/m^3
    density=1000; %in kg/m^3
    cost_per_kwh=.08;
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
        cost_per_kwh = .05;
    end
    pump_cost_perunit=density * 9.81 * (((f*(s1+aquiferDepth)^2*v^2)/(2*9.81*D))+(s1+aquiferDepth)) * conversion_factor ...
        * (100/percent_eff) * cost_per_kwh;

    % Expansion costs
    switch a2
        case 0
            expansionCost = 0;
        case 1
            capex = (-0.0014 * water.desal_capacity_expansion.small/365 + 1362.5) * water.desal_capacity_expansion.small/365;
            expansionCost = capex * discountFactor;
        case 2
            capex = (-0.0014 * water.desal_capacity_expansion.large/365 + 1362.5) * water.desal_capacity_expansion.large/365;
            expansionCost = capex * discountFactor;
    end

% Third: Calculate supply

    % Always use groundwater if we have it 
    minjur_supply = minjur_capacity;
    
    % Calculate water supplied by desal exp
    % exp_supply = max(0, demand - minjur_supply - water.desal_capacity_initial - othergw_supply);
    exp_supply = max(0, demand - minjur_capacity);
    exp_supply = min(exp_supply, exp_capacity);
    
    % Calculate desal marginal costs
%     opexSmall = (-0.0001 * water.desal_capacity_expansion.small/365 + 186.76) * water.desal_capacity_expansion.small/365;
%     opexLarge = (-0.0001 * water.desal_capacity_expansion.large/365 + 186.76) * water.desal_capacity_expansion.large/365;
    opex = 170 / 365;  % $/m^3/y

    % Calculate shortage
    shortage = max(0, demand - capacity);
   

% Costs include shortage costs, expansion costs, and pumping costs
    shortageCost = shortage * costParam.shortage_cost * discountFactor;
    pumpingCost =  pump_cost_perunit * a1 * minjur_supply * discountFactor;
    marginalDesalCost = exp_supply * costParam.marginal_cost;
    cost = shortageCost + expansionCost + pumpingCost + marginalDesalCost;


end

