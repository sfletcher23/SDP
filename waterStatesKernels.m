function [s_gw, gw_M, kernel, index_T_S_samples] = waterStatesKernels(gwParam, s_pop, water)

N = 2;
%% State and Action Definitions for Groundwater 


% Import groundwater well data and aquifer properties
% T in units m^2/sec, drawdown and depth in unit meters
groundwaterWells = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A6:I19', 'ReadVariableNames', true, 'ReadRowNames', true);
aquifer = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A25:J27', 'ReadVariableNames', true, 'ReadRowNames', true);

% For now, just model first well with input demand per capita
a = 1;  % Selects first well
aquiferNames = aquifer.Properties.RowNames;
pumpLocation = [groundwaterWells.Lattitude(a) groundwaterWells.Longitude(a)];
observeLocation = 10;

% Get min and max demand
minDemand = demand(water, s_pop(1));
maxDemand = demand(water, s_pop(end));

% Theis Paramters
pumpStep = 0;
QTheisMax = maxDemand;
QTheisMin = zeros(N,1);
QTheisMin(1) = minDemand;
locUnits = 'meters';

%% Min and Max drawdown
    
Tmax = aquifer.TransmissivityUpperBound(aquiferNames{a}) * 60 * 60 *24 * 365; % [ m^2/year]
Tmin = aquifer.TransmissivityLowerBound(aquiferNames{a}) * 60 * 60 *24 * 365; % [ m^2/year]

% Deterministic storativity
S = aquifer.Storativity(aquiferNames{a});   

% Using Theis restricted didn't work because of only 1 pumpstep. In
% general, need to figure out how to deal with 

outputmax = theis(QTheisMax, pumpStep, Tmin, S, pumpLocation, observeLocation, locUnits, 1,1);
drawdownMaxAnnual = max(outputmax);
outputmin = theis(QTheisMin, pumpStep, Tmax, S, pumpLocation, observeLocation, locUnits, 1,1);
drawdownMinAnnual = min(outputmin);


%% State and Action Definitions


%  Calculate discretization size
    % avgDiff = mean([diff(outputmin) diff(outputmax)]);  % Average difference between drawdown impacts
    % step = ceil(avgDiff);
    %step = round(ceil(drawdownMinAnnual)*1.5);
    step = ceil(drawdownMinAnnual);

% Define states: 
s_gw =[ 0: step: gwParam.depthLimit]';
gw_M = length(s_gw); 

% Round up annual max
drawdownMaxAnnual = (floor(drawdownMaxAnnual / step) + 1) * step ;

%% Calculate kernel functions
% Calculate kernel functions measuring the response to a unit impulse of
% pumping from each well over time. 
% Kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix
% T_S_pairs is a vector of samples from T and S 

[kernel, T_S_pairs] = gen_kernel(groundwaterWells, aquifer, N);
% [numSamples, ~] = size(T_S_pairs);
% kernel = reshape(kernel(:,2,:), [1,numSamples]);

    % Test that kernel functions provide reasonable drawdowns for demand
    testmin = min(max(kernel * minDemand));
    testmax = max(max(kernel * maxDemand));
    test1 = testmin > gwParam.stepSize;
    test2 = testmax < max(s_gw);
    if ~test1
        error('Kernels generate unrealistically low drawdown')
    end
    if ~test2
        error('Kernels generate unrealistically high drawdown')
    end


%% Calculate pt(k), the distribution over parameter vector k in time t
    % and generate samples of k for each t

index_T_S_samples = zeros(gwParam.sampleSize,N);
for t = 1:N
    % Define pt(k)
        [T_S_pair_cdf] = gen_param_dist(T_S_pairs, t);
    % Generate samples from pt(K)
        p = rand(1,gwParam.sampleSize);
        index_T_S_samples(:,t) = arrayfun(@(x) find(x < T_S_pair_cdf,1), p);
        clear T_S_pair_cdf p
end



end