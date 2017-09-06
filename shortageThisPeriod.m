function [shortage, supply, demand, gw_supply, exp_supply] = shortageThisPeriod(a1, s1, s2, water, demand, s_gw, gwParam)
% Calculates shortage in current period gives states, actions, and demand

% Groundwater supplied this period
if a1 == 0 || s1 == 200 || s1 == -1
    gw_supply = 0;
else
    gw_supply = gwParam.pumpingRate;
end
	

% Exp desal supply
exp_supply = s2; % assume use 100% of desal capacity

% Total water supplied is sum of groundwater and desalinated water
supply = gw_supply + water.desal_capacity_initial + exp_supply; 
 
% Calculate shortage
shortage = max(0, demand - supply);


end

% update this for gw supply to reflect withdawal demand insted of a1
 
