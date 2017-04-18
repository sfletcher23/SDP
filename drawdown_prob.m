function [dd_prob dd_values] = drawdown_prob(gwParam, demand_range, groundwaterWells, aquifer, drawdownMaxAnnual)
%% Theis Parameters

a = 1;  % first well for now
pumpLocation = [groundwaterWells.Lattitude(a) groundwaterWells.Longitude(a)];
observeLocation = 10;
pumpStep = 0;
time = 1;
locUnits = 'meters';

%% Simulate transmissivity and aquifer depth

aquiferNames = aquifer.Properties.RowNames;
distT = makedist('triangular', aquifer.TransmissivityLowerBound(aquiferNames{a}), ...
    aquifer.TransmissivityAverage(aquiferNames{a}), aquifer.TransmissivityUpperBound(aquiferNames{a}));
randT = random(distT,gwParam.sampleSize,1);
T = randT * 60 * 60 *24 * 365; % [ m^2/year]
% distAquiferDepth = makedist('triangular', aquifer.DepthHigh(aquiferNames{a}), ...
%     aquifer.DepthAverage(aquiferNames{a}), aquifer.DepthLow(aquiferNames{a}));
% randDepth = random(distAquiferDepth,gwParam.sampleSize,1);

% Deterministic storativity
S = ones(gwParam.sampleSize,1)*aquifer.Storativity(aquiferNames{a});   


%% Get drawdown for each transmissivity and aquifer depth sample

drawdown = zeros(gwParam.sampleSize, length(demand_range));
for j = 1:length(demand_range)
     drawdown(:,j) = theis( demand_range(j), pumpStep, T, S, pumpLocation, observeLocation, locUnits, time, gwParam.sampleSize );
end

dd_values = 0: gwParam.stepSize : drawdownMaxAnnual;
roundedDrawdown = round2x(drawdown, dd_values);

dd_prob = zeros(length(demand_range),length(dd_values));
for i = 1:length(demand_range)
    dd_prob(i,:) = histcounts(roundedDrawdown(:,i), [-5 dd_values], 'Normalization', 'probability');
end

% Test that all rows sum to 1
err = sum(dd_prob,2) - 1;
margin =  1E-4;
if err > margin
    error('Invalid drawdown distribution')
end


% How to handle time. Add in state space? Time in Theis measured as time since
% start of pumping. Assume that starting is in equilibrium? Or not? 
% Currently: just adding drawdown to previous drawdown after 1 period

end

