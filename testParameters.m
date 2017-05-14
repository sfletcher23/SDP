 %% Test paramters
 
 %
 
 % Initial, max, and min annual demand
 demand_initial = demand(water, popParam.pop_initial, water.demandFraction);
 demand_min = min(demand_range);
 demand_max = max(demand_range);
 
 
 % Initial annual capacity
 gwCapacity = gwParam.pumpingRate;
 desalCapacity = water.desal_capacity_initial;
 expansionCapacity = water.desal_capacity_expansion;
 

% Plot costs 
% min costs
shortage_min_gw_exp = max(demand_min - (gwCapacity + desalCapacity + expansionCapacity), 0);
shortage_min_gw_noexp = max(demand_min - (gwCapacity + desalCapacity), 0);
shortage_min_nogw_exp = max(demand_min - (desalCapacity + expansionCapacity), 0);
shortage_min_nogw_noexp = max(demand_min - (desalCapacity), 0);

shortage_max_gw_exp = max(demand_max - (gwCapacity + desalCapacity + expansionCapacity), 0);
shortage_max_gw_noexp = max(demand_max - (gwCapacity + desalCapacity), 0);
shortage_max_nogw_exp = max(demand_max - (desalCapacity + expansionCapacity), 0);
shortage_max_nogw_noexp = max(demand_max - (desalCapacity), 0);

cost_min_gw_exp = [shortage_min_gw_exp * costParam.shortage_cost* N, costParam.expansion_cost, costParam.pumping_cost * N];
cost_min_gw_noexp = [shortage_min_gw_noexp * costParam.shortage_cost* N, 0 , costParam.pumping_cost * N];
cost_min_nogw_exp = [shortage_min_nogw_exp * costParam.shortage_cost* N , costParam.expansion_cost, 0];
cost_min_nogw_noexp = [shortage_min_nogw_noexp * costParam.shortage_cost* N, 0 0];

cost_max_gw_exp = [shortage_max_gw_exp * costParam.shortage_cost* N, costParam.expansion_cost, costParam.pumping_cost * N];
cost_max_gw_noexp = [shortage_max_gw_noexp * costParam.shortage_cost* N, 0, costParam.pumping_cost * N];
cost_max_nogw_exp = [shortage_max_nogw_exp * costParam.shortage_cost* N, costParam.expansion_cost, 0];
cost_max_nogw_noexp = [shortage_max_nogw_noexp * costParam.shortage_cost* N, 0, 0];

labels = {'Pumping, expanion', 'Pumping, no expansion', 'No pumping, expansion', 'No pumping, no expansion'};


 % Plot demand vs. capacity
 figure;
 y = [gwCapacity desalCapacity expansionCapacity; 0 0 0];
bar(y,'stacked')
legend('GW', 'desal intial', 'expansion')
hold on
scatter([2 2 2], [demand_initial, demand_max, demand_min]);


% plot costs
figure
subplot(2,1,1)
b = bar([cost_min_gw_exp; cost_min_gw_noexp; cost_min_nogw_exp; cost_min_nogw_noexp],'stacked');
set(gca,'XTick', 1:4, 'XTickLabel', labels);
title('min demand')
legend('shortage costs', 'exp costs', 'pumping costs')
subplot(2,1,2)
b = bar([cost_max_gw_exp; cost_max_gw_noexp; cost_max_nogw_exp; cost_max_nogw_noexp],'stacked');
set(gca,'XTick', 1:4, 'XTickLabel', labels);
title('max demand')
legend('shortage costs', 'exp costs', 'pumping costs')
