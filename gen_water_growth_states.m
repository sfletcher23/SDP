function [s_gw, gw_M, drawdownMaxAnnual, step] = gen_water_growth_states(gwParam, N, demand_range, water, groundwaterWells, aquifer)
% Takes input paramters and length of model and generates groundwater state space 


% For now, just model first well with input demand per capita
a = 1;
aquiferNames = aquifer.Properties.RowNames;
pumpLocation = [groundwaterWells.Lattitude(a) groundwaterWells.Longitude(a)];
observeLocation = 10;

% Theis Paramters
pumpStep = 0;
QTheisMax = max(demand_range);
QTheisMin = zeros(N,1);
QTheisMin(1) = min(demand_range);
locUnits = 'meters';

%% Min and Max drawdown
    
Tmax = aquifer.TransmissivityUpperBound(aquiferNames{a}) * 60 * 60 *24 * 365; % [ m^2/year]
Tmin = aquifer.TransmissivityLowerBound(aquiferNames{a}) * 60 * 60 *24 * 365; % [ m^2/year]

% Deterministic storativity
S = aquifer.Storativity(aquiferNames{a});   

% Using Theis restricted didn't work because of only 1 pumpstep. In
% general, need to figure out how to deal with 

output = theis(QTheisMax, pumpStep, Tmin, S, pumpLocation, observeLocation, locUnits, 1:N,1);
drawdownMaxAnnual = max(output);
output = theis(QTheisMin, pumpStep, Tmax, S, pumpLocation, observeLocation, locUnits, 1:N,1);
drawdownMinAnnual = min(output);

drawdownMax = drawdownMaxAnnual * N;
drawdownMin = drawdownMinAnnual * N;

%% State and Action Definitions

% Which is limiting: pumping or aquifer depth?
limit = min(drawdownMax, gwParam.depthLimit);

%  Use min annual drawdown as discretization size
step = floor(drawdownMinAnnual);

% Define states: 
s_gw = 0: step: limit;
gw_M = length(s_gw); 

% Round up annual max
drawdownMaxAnnual = (floor(drawdownMaxAnnual / step) + 1) * step ;


end
