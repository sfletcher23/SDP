%% Test kernel function method
% Calculate drawdown for a specific pumping scenario using both methods and
% make sure they are the same

%% Import groundwater well data and aquifer properties
% T in units m^2/sec, drawdown and depth in unit meters

groundwaterWells = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A6:I19', 'ReadVariableNames', true, 'ReadRowNames', true);
aquifer = readtable('Water Decision Model Data Inputs.xlsx', 'Sheet', 4, ...
    'Range', 'A25:J27', 'ReadVariableNames', true, 'ReadRowNames', true);

%% Define ppumping inputs

a = 1:5;  % well rows
pumpLocation = [groundwaterWells.Lattitude(a) groundwaterWells.Longitude(a)];
observeLocation = 10;
t = 1:1:10;
pumpStep = t;
locUnits = 'meters';
Q = repmat(groundwaterWells.WithdrawalsIn2010(a),1,length(t))'*1E6;

%% Define parameter space

T_length = 10;
S_length = 5;
numParameterValues = T_length*S_length;

aquiferNames = aquifer.Properties.RowNames;
T_lower = aquifer.TransmissivityLowerBound('Minjur')* 60 * 60 *24 * 365; % [ m^2/year]
T_upper = aquifer.TransmissivityUpperBound('Minjur')* 60 * 60 *24 * 365; % [ m^2/year]
S_lower = aquifer.Storativity('Minjur') / 10;
S_upper = aquifer.Storativity('Minjur') * 10;
T_temp = linspace(T_lower,T_upper,T_length);
S_temp = linspace(S_lower, S_upper,S_length);

Stemp = repmat(S_temp,T_length,1);
[row, col] = size(Stemp);
Stemp2 = reshape(Stemp,row*col,1);
Ttemp = repmat(T_temp,1,S_length)';
T_S_pairs = [Ttemp Stemp2];
[runs,~] = size(T_S_pairs);
T = T_S_pairs(:,1);
S = T_S_pairs(:,2);

%% Calculate drawdown from multiple pumping wells over time using Theis
% Drawdown is a [numObserve x numTime x runs] matrix

[ drawdownTheis ] = theis(Q, pumpStep, T, S, pumpLocation, observeLocation, locUnits, t, runs );


%% Calculate drawdown from multiple pumping wells over time using kernels
% Each kernel function measures drawdown at a set of observation locations
% for unit pumping impulse from a single pumping well
% Kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix

[numPumpWells, ~] = size(pumpLocation);
numObserve = numPumpWells;
numTime = length(t);
kernel = zeros(numObserve,numTime,numParameterValues,numPumpWells);
pumpStep = t;
for i = 1:numPumpWells
    Qkernel = zeros(numTime,numPumpWells);
    Qkernel(1,i) = 1; 
    drawdown = theis( Qkernel, pumpStep, T, S, pumpLocation, observeLocation, locUnits, t, numParameterValues );
    kernel(:,:,:,i) = drawdown; 
end

drawdownKernel = kernel2drawdown(Q, kernel, t);

%% Compare drawdownTheis with drawdownKernel

difference = abs(drawdownTheis - drawdownKernel);
tolerance = 1e-3;
indexError = difference > tolerance;
sumError = sum(sum(sum(indexError)));
if sumError > 0 
    error('Methods are not equal')
elseif sumError == 0
    disp('Test passed')
end



