function [s_gw, gw_M, drawdownMaxAnnual, step] = gen_water_growth_states(gwParam, N, demand_range, water, groundwaterWells, aquifer)

% Takes input paramters and length of model and generates groundwater state space 

% Inputs
    % gwParam is a struct with input paramters on groundwater
    % N is the time horizon of model
    % demand_range is a vector of possible groundwater demand values based on population discretization
    % water is a struct with input paramters on desal and other infrastructure
    % groundwater wells is a table of data on the groundwater wells in the two major aquifers
    % aquifer wells is a table of data on the major aquifers
    
% Outputs
    % s_gw is a vector with the possible drawdown levels in the groundwater well
    % gw_M is the number of states in s_gw
    % drawdownMaxAnnaul is the largest possible drawdown in a single year
    % step is the discretization size of the

% For now, just model first well with input demand per capita
a = 1;  % Selects first well
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

outputmax = theis(QTheisMax, pumpStep, Tmin, S, pumpLocation, observeLocation, locUnits, 1:N,1);
drawdownMaxAnnual = max(outputmax);
outputmin = theis(QTheisMin, pumpStep, Tmax, S, pumpLocation, observeLocation, locUnits, 1:N,1);
drawdownMinAnnual = min(outputmin);

drawdownMax = drawdownMaxAnnual * N;
drawdownMin = drawdownMinAnnual * N;

%% State and Action Definitions

% Which is the upper limit on drawodnw: pumping or aquifer depth?
limit = min(drawdownMax, gwParam.depthLimit);

%  Calculate discretization size
    % avgDiff = mean([diff(outputmin) diff(outputmax)]);  % Average difference between drawdown impacts
    % step = ceil(avgDiff);
step = ceil(drawdownMinAnnual / 2);

% Define states: 
s_gw = 0: step: limit;
gw_M = length(s_gw); 

% Round up annual max
drawdownMaxAnnual = (floor(drawdownMaxAnnual / step) + 1) * step ;


end
