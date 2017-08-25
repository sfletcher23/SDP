function [s_gw, gw_M] = gen_water_growth_states(gwParam, N, demand_range, water, groundwaterWells, aquifer)

% Takes input paramters generates groundwater state space 
limit = gwParam.depthLimit;

%  Calculate discretization size
step = 1;

% Define states: 
s_gw = 0: step: limit;
gw_M = length(s_gw); 

end
