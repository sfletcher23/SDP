function [ cost, shortageCost, expansionCost, pumpingCost, marginalDesalCost, shortage, capacity, minjur_supply, exp_supply] ...
    = supplyAndCost_old( a1, a2, s1, s2, costParam, water, gwParam, t, demand, capacityDelay, exp_vectors, makePlot)


% Calculate the capacity, supply, and cost in the current period given
% acitons and state

% First: Calculate capacity available

    % Groundwater capacity this period
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
    density=1000; %in kg/m^3
    cost_per_kwh=.08;
    conversion_factor=2.777778e-7;
    percent_eff=80;
    f=0.013;
    D=0.508;
    v=0.24669;
    aquiferDepth = 1200;
    startingDepth = gwParam.startingHead;
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
            capex = (-189.1 * log(water.desal_capacity_expansion.small/365) + 3385) * water.desal_capacity_expansion.small/365;
            expansionCost = capex * discountFactor;
        case 2
            capex = (-189.1 * log(water.desal_capacity_expansion.large/365) + 3385) * water.desal_capacity_expansion.large/365;
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
    
    
    if makePlot
        figure;
        f_pump = f;
        D_pump = D;
        v_pump = v;
        density_pump = density;
        percent_eff_pump = percent_eff;
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
        
        fig = figure; 
        colors = get(fig, 'defaultAxesColorOrder');
        set(fig,'defaultAxesColorOrder',[colors(1,:); colors(2,:)]);
        yyaxis left
        bl = bar([costParam.shortage_cost 0; pc(1) pc(end) - pc(1); 0 0; costParam.marginal_cost 0; 0 0], 'stacked');
        ylim([0 10])
        ylabel('Marginal Costs [$/cm]')
        yyaxis right
        br = bar([0; 0 ;0 ;0 ; expansionCost/1E6], 'FaceColor', colors(2,:));
        ylim([0 700])
        ylabel('Capex [Million $]')
        set(gca, 'XTickLabel', {'shortage', 'pumping', 'brackish', 'desal opex', 'desal capex'})
        bl(1).FaceColor = colors(1,:);
        bl(2).FaceColor = colors(1,:);
        bl(2).FaceAlpha = 0.4;
        b = findall(fig,'Type', 'Bar');
        set(b, 'BarWidth', 0.6)
        hold on
        line([4.5 4.5], [0 700], 'Color', 'k')
        title('Cost Assumptions')
        expansionCost/1E6
    end


end