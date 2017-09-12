function [shortage, capacity, demand, minjur_supply, exp_supply, othergw_supply] = shortageThisPeriod(a1, s1, s2, water, demand, s_gw, gwParam)
% Calculates shortage in current period gives states, actions, and demand

% Groundwater supplied this period
othergw_supply = gwParam.otherPumpingRate;
if a1 == 0 || s1 == 200 || s1 == -1
    minjur_supply = 0;
else
    minjur_supply = gwParam.pumpingRate;
end
	

% Exp desal supply
exp_capacity = s2; % assume use 100% of desal capacity

% Total water capacity is sum of groundwater and desalinated water
capacity = minjur_supply + water.desal_capacity_initial + exp_capacity + othergw_supply  ; 
 
% Calculate shortage
shortage = max(0, demand - capacity);

% Calculate water supplied by desal exp
exp_supply = max(0, demand - minjur_supply - water.desal_capacity_initial);
exp_supply = min(exp_supply, exp_capacity);

end


 
