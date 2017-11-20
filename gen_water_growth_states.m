function [s_gw, gw_M] = gen_water_growth_states(gwParam)

% Takes input paramters generates groundwater state space 
limit = 588 + 337;

% If depth limit, cap state space
if gwParam.depthLimit
    limit = gwParam.depthLimit;
end

%  Calculate discretization size
step = 0.5;

% Define states: 
s_gw = 0: step: limit;
s_gw = [-1 s_gw]; % This is absorbing state where can't pump anymore
gw_M = length(s_gw); 

end
