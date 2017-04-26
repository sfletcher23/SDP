function [ drawdown ] = theis( Q, pumpStep, T, S, pumpLocation, observeLocation, locUnits, t, runs )
% Input location and parameters for each pumping well and time 
% Output drawndon in each observation well over time 

    % Assumes: Confined aquifer, no recharge, observations wells sufficiently 
    % far away from pumping wells
    
    % Drawdown is a [numObserve x numTime x runs] matrix
    % that gives the drawdown in each obseration well over time.

% Throw errors if incorrect input sizes

% Q is a [time X numPumpWells] matrix


% Number parameters
[numPumpWells,~] = size(pumpLocation);  % Number of pumping wells
[numObserve, ~] = size(observeLocation);    % Number of observations
numTime = length(t);    % Number of time steps
numPumpStep = length(pumpStep); 

% Input for observeLocation can be either a set of locations (coordinates)
% or a single radial distance from well. If a single radial distance r (assumed
% in meters), find a point r distance away from each (in arbitrary direction). 
if  size(observeLocation) == [1 1]   
    if strcmp(locUnits,'meters') 
        r = observeLocation;
        observeLocation = pumpLocation + [ones(numPumpWells,1)*r  zeros(numPumpWells,1) ]; 
        [numObserve, ~] = size(observeLocation);
    elseif strcmp(locUnits, 'degrees')
        r = distdim(observeLocation,'meters','degrees');
        observeLocation = pumpLocation + [zeros(numPumpWells,1) ones(numPumpWells,1)*r   ]; 
        [numObserve, ~] = size(observeLocation);
    else
        error('Incorrect locUnits for Theis: must be meters or degrees')
    end
end
  
%% Calculate distance from pumping to observation
% R is a [numObserve x numPumpWells] symmetric matrix defining the distance
% between each obsevation location and each pumping well. If locations
% input as degrees instead of meters, convert distances to meters.

R = zeros(numObserve,numPumpWells);

if strcmp(locUnits,'meters') 
    for i = 1:numObserve
        for j = 1:numPumpWells
            R(i,j) = sqrt( (pumpLocation(j,1) - observeLocation(i,1)) .^2 + ...
                (pumpLocation(j,2) - observeLocation(i,2)) .^2 );
        end
    end
    
elseif strcmp(locUnits, 'degrees')
    for i = 1:numObserve
        for j = 1:numPumpWells
            R(i,j) = distdim(distance(pumpLocation(j,:), observeLocation(i,:)),'degrees', 'meters');
        end
    end
    
else
    error('Incorrect locUnits for Theis: must be meters or degrees')
end

    % RepR2 is a [numObserve x numPumpWells x numTime x runs] matrix that gives
    % the distance squared between each observation location and each pumping
    % well and replicates it for each time step for future calculations
    repR2 = repmat(R.^2, [1,1,numTime,runs]);

%% Calculate drawdown using Theis
% Repeat drawdown calculation for every pumping time step

deltaQ = Q(1,:);
drawdownStep = zeros([numObserve, numTime, numPumpStep, runs]);

for k = 1:numPumpStep
    
    % Calculate tk - set times before current step equal to zero so new
    % image well doesn't impact them.  
    tk = t - pumpStep(k);
    index = tk < 0;
    tk(index) = 0;

    % Calculate well function for superposed well. W is a [numObserve x 
    % numPumpWells x numTime x runs] matrix, which describes the well function
    % for each observation location from each well in every time step.
    rep_t = repmat(tk', [1,numObserve,numPumpWells,runs]);  % [numTime x numObserve x numPumpWells x runs]
    rep_t = permute(rep_t, [2 3 1 4]);  % [numObserve x numPumpWells x numTime x runs]
    rep_T = repmat(T', [numObserve,1,numPumpWells,numTime]);    % [numObserve x runs x numPumpWells x numTime]
    rep_T = permute(rep_T, [1 3 4 2]);  % [numObserve x numPumpWells x numTime x runs]
    rep_S = repmat(S', [numObserve,1, numPumpWells, numTime]);   % [numObserve x runs x numPumpWells x numTime]
    rep_S = permute(rep_S, [1 3 4 2]);  % [numObserve x numPumpWells x numTime x runs]
    U = repR2 .* rep_S ./ (4 * rep_T .* rep_t);
    W = arrayfun(@(x) expint(x), U);
    W(isnan(W)) = 0;
    
    % Update deltaQ, the difference in pumping rate since previous
    % step. [1 X numPumpWells]
    % RepQ is a [numObserve x numPumpWells x numTime X runs] matrix that
    % gives deltaQ from each pumping well for each observation well at
    % all times t.
    if k > 1
        deltaQ = Q(k,:) - Q(k-1,:);
    end
    repQ = repmat(deltaQ,[numObserve, 1 , numTime, runs]);

    % Calculate drawdown contribution in each observation well at time 
    % t, using repQ and well function. Sum drawdown contributions over
    % all pumping wells 
    % Drawdown is a [numObserve x numTime x runs] matrix
    % that gives the drawdown in each obseration well over time.
    drawdownContributions = W .* repQ ./ (4*pi*rep_T); % [numObserve x numPumpWells x numTime X runs]
    drawdownTemp = sum(drawdownContributions,2); % sum over pumpnig wells
    drawdownStep(:,:,k,:) = reshape(drawdownTemp, [numObserve, numTime, runs]);
end
drawdown = reshape(sum(drawdownStep,3),[numObserve numTime runs]);

