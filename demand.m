function [ demand_range ] = demand( water, pop, t)

% Caluclates water emand, given water demand per capita, population (vector or
% scalar), and the fraction of total demand 

    % If s_pop is a vector, outputs vector wiht possible demand levels. If
    % s_pop is a scalar, outputs scalar wtih demand.

demand_range = (water.demandPerCapita(t) * 365/1000)  ...    % m^3/p/year
    .* pop * 1E6 ...   % p
    * water.demandFraction;


end