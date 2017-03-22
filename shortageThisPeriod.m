function [shortage, supply] = shortageThisPeriod(a1, a2, s1, s2, water)
% Calculates shortage in current period gives states, actions, and demand

% Groundwater supplied this period given by a1
gw_supply = a1;

% Total water supplied is sum of groundwater and desalinated water
supply = gw_supply + water.desal_capacity_initial + ...
    water.desal_capacity_expansion * (s2-1); % desal if available

% Calculate shortage
shortage = max(0, water.demand - supply);


end

