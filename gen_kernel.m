function [kernel, T_S_pairs] = gen_kernel(groundwaterWells, aquifer, N)

% Each kernel function measures drawdown at a set of observation locations
% for unit pumping impulse from a single pumping well

% Kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix

%% Theis Parameters

a = [1 ];  % well rows
pumpLocation = [groundwaterWells.Lattitude(a) groundwaterWells.Longitude(a)];
observeLocation = 10;
pumpStep = 0;
time = 1:1:N;
locUnits = 'meters';


%% Discretize parameter space

T_length = 100;
S_length = 50;
numParameterValues = T_length*S_length;

aquiferNames = aquifer.Properties.RowNames;
T_lower = aquifer.TransmissivityLowerBound('Minjur')* 60 * 60 *24 * 365; % [ m^2/year]
T_upper = aquifer.TransmissivityUpperBound('Minjur')* 60 * 60 *24 * 365; % [ m^2/year]
S_lower = aquifer.Storativity('Minjur') / 10;
S_upper = aquifer.Storativity('Minjur') * 10;
T = linspace(T_lower,T_upper,T_length);
S = linspace(S_lower, S_upper,S_length);

Stemp = repmat(S,T_length,1);
[row, col] = size(Stemp);
Stemp2 = reshape(Stemp,row*col,1);
Ttemp = repmat(T,1,S_length)';
T_S_pairs = [Ttemp Stemp2];


%% Run Model for every parameter (T,S pair) and calculate kernel value

% Drawdown is a [numObserve x numTime x numParameterValues] matrix
    % that gives the drawdown in each obseration well over time.
    
% In each Theis run, only 1 of the pumping wells is on, so loop over all
% the pumping wellls

% Each kernel function measures drawdown at a set of observation locations
% for unit pumping impulse from a single pumping well
% Kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix


[numPumpWells, ~] = size(pumpLocation);
numObserve = numPumpWells;
numTime = length(time);
kernel = zeros(numObserve,numTime,numParameterValues,numPumpWells);
pumpStep = time(1:2);
for i = 1:numPumpWells
    Q = zeros(numTime,numPumpWells);
    Q(1,i) = 1; 
    drawdown = theis( Q, pumpStep, T_S_pairs(:,1), T_S_pairs(:,2), pumpLocation, observeLocation, locUnits, time, numParameterValues );
    kernel(:,:,:,i) = drawdown; 
end


end
