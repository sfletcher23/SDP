function [ cost, shortageCost, expansionCost, pumpingCost, marginalDesalCost, shortage, capacity, minjur_supply, exp_supply, othergw_supply ] ...
    = supplyAndCost( a1, a2, s1, s2, costParam, water, gwParam, t, demand, capacityDelay, exp_vectors)


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
    if capacityDelay
        s2index = linIndex2VecIndex(s2, exp_vectors);
        exp_capacity = exp_vectors{1}(s2index(1));
    end

    % Total water capacity is sum of groundwater and desalinated water
    % capacity = minjur_supply + water.desal_capacity_initial + exp_capacity + othergw_supply;
    capacity = minjur_capacity + exp_capacity ; 

% Second: Calculate expansion costs and marginal costs

    % Discount factor
    discountFactor = 1/((1+costParam.discount_rate)^t);

    % Pumping cost in $/m^3
    density_pump=1010; %in kg/m^3
    cost_per_kwh=.10; %$/kWh
    conversion_factor=2.777778e-7;  % kWh/J
    percent_eff_pump=80;
    f_pump=0.04;
    D_pump=0.2;
    v_pump=18;
    startingDepth = 300;
    brackish_cost_min = .3;
    brackish_cost_max = .6; % At 300 m drawdown
    brackish_cost =  brackish_cost_min + (brackish_cost_max - brackish_cost_min) * s1/300;
    effective_head_pump = (s1 + startingDepth) + f_pump * (s1+startingDepth)/D_pump * v_pump^2 / (2*9.81);
    pump_cost_perunit=density_pump * 9.81 * effective_head_pump * conversion_factor ...
        * (100/percent_eff_pump) * cost_per_kwh;

    
    % Pipeline costs
    percent_eff_pipe=90;
    density_pipe=1000;
    f_pipe=0.01;
    D_pipe=1.6;
    v_pipe=4;
    L = 466000;
    delta_h = 615;
    effective_head_pipe = delta_h + f_pipe * L/D_pipe * v_pipe^2 / (2*9.81);
    pipe_cost_perunit = density_pipe * 9.81 * effective_head_pipe * conversion_factor ...
        * 100/percent_eff_pipe * cost_per_kwh;
    
    
    if false
        figure;
        for s1 = 1:400
            Hf(s1) = (s1 + startingDepth) + f_pump * (s1+startingDepth)/D_pump * v_pump^2 / (2*9.81);
            pc(s1) = density_pump * 9.81 * Hf(s1) * conversion_factor ...
            * (100/percent_eff_pump) * cost_per_kwh;
        end
        subplot(1,2,1)
        plot((1:400)+ startingDepth, pc)
        xlabel('drawdown [m]')
        ylabel('pumping cost [$/m3]')
        title('Pumping costs')
        subplot(1,2,2)
        plot((1:400) + startingDepth, Hf)
        xlabel('drawdown [m]')
        ylabel('effective depth [m]')
        title('Effective pumping depth')
    end

    % Expansion costs
    expansionCost = [];
    switch a2
        case 0
            expansionCost = 0;
        case 1
            capex = (-189.1 * log(water.desal_capacity_expansion.small/365) + 3385) * water.desal_capacity_expansion.small/365;
            expansionCost = capex * discountFactor;
%             expansionCost = expansionCost * 1.2; % Just to test
        case 2
            capex = (-189.1 * log(water.desal_capacity_expansion.large/365) + 3385) * water.desal_capacity_expansion.large/365;
            expansionCost = capex * discountFactor;
%             expansionCost = expansionCost * 1.2; % Just to test
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
    pumpingCost =  (pump_cost_perunit + brackish_cost) * a1 * minjur_supply * discountFactor;
    marginalDesalCost = exp_supply * costParam.marginal_cost;
    cost = shortageCost + expansionCost + pumpingCost + marginalDesalCost;


end

