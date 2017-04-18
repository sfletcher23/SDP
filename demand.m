function [ demand_range ] = demand( water, s_pop, fraction)

% Caluclates demand, given water demand per capita, population (vector or
% scalar), and the fraction of total demand 

demand_range = water.demandPerCapita * 365/1000  ...    % m^3/p/year
    * s_pop * 1E6 ...   % p
    * fraction;


end