function [ drawdown, QRestricted, unmetDemand] = theisRestricted( QTheis, T, S, pumpLocation, observeLocation, time, locUnits, drawdownLimit, runs )
% Input location and parameters for each pumping well and time 
% Output drawndon in each observation well over time using Theis equation

% Time units: can be anything as long as pumping rate and time are the same
% unit. Check on this. 

% Linear superposition is used to represent the impacts of multiple wells
% on each observation location.
% Linear superposition is also used to represent changes in pumping rate
% over time by adding a well with rate deltaQ every time rate changes.

    % Assumes: Confined aquifer, no recharge, observations wells sufficiently 
    % far away from pumping wells
    
    % If observeLocation is a matrix, represents coordinates of observation
    % wells. If scalar, represents distance away from pumping wells.
   
% % Throw errors if incorrect input sizes
%     % Q should be same length as time
%     if length(QTheis) ~= length(time)
%         error('QTheis must be same length as time')
%     end

%% Parameters 
% Number parameters
[numPumpWells, ~] = size(pumpLocation);  % Number of pumping wells
[numObserve, ~] = size(observeLocation);    % Number of observations
numTime = length(time);    % Number of time steps
numPumpStep = length(time); % Number of times rate of pumping changes
QInput = QTheis;

[row, col] = size(QTheis);    
QRestricted = zeros(row,col,runs);
unmetDemand = zeros(row,col,runs);
drawdown = zeros(numPumpWells,numTime,runs); 
    
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

    % RepR2 is a [numObserve x numPumpWells x numTime] matrix that gives
    % the distance squared between each observation location and each pumping
    % well and replicates it for each time step for future calculations
    repR2 = repmat(R.^2, [1,1,numTime]);


%% Calculate drawdown in well for each time step

for n = 1:runs
    
    % DeltaQ is a [numPumpWells] length array that gives the change in 
    % pumping rate since previous change, to be used for superposition 
    % drawdownStep gives
    QTheisRun = QTheis;
    deltaQ = QTheis(1,:);
    drawdownStep = zeros([numObserve,numTime,numPumpStep]);

    % For each change in pumping rate, add a superposed set of well and 
    % calculate drawdown contribution in each observation well over time
    for t = 1:length(time)

        % Calculate tk - set times before current step equal to zero so new
        % image well doesn't impact them.
        tk = time - t;
        indexTime = tk < 0;
        tk(indexTime) = 0;

        % Calculate well function for superposed well. W is a [numObserve x 
        % numPumpWells x numTime] matrix, which describes the well function
        % for each observation location from each well in every time step.
        
        
        t3 = zeros(1,1,length(tk));
        t3(1,1,:) = tk;
        repT = repmat(t3, [numObserve,numPumpWells,1]);
        U = repR2* S(n) ./ (4 * T(n) * repT);
        W = arrayfun(@(x) expint(x), U);
        W(isnan(W)) = 0;

        % Update deltaQ, the difference in pumping rate since previous
        % step. RepQ is a [numObserve x numPumpWells x numTime] matrix that
        % gives deltaQ from each pumping well for each observation well at
        % all times t.
     
        if t > 1
            deltaQ = QTheisRun(t,:) - QTheisRun(t-1,:);
        end
        repQ = repmat(deltaQ,[numObserve, 1 , numTime]);

        % Calculate drawdown contribution in each observation well at time 
        % t, using repQ and well function. Sum drawdown contributions over
        % all pumping wells. drawdown is a [numObserve x numTime] matrix
        % that gives the drawdown in each obseration well over time.
        drawdownContributions = W .* repQ / (4*pi*T(n));
        drawdownTemp = sum(drawdownContributions,2);
        drawdownStep(:,:,t) = reshape(drawdownTemp, [numObserve, numTime]);
        drawdownRun = sum(drawdownStep,3);
        

        % If drawdown exceeds drawdownThreshold, reduce pumping by factor of 1.2. 
        % Update QTheisRun (pumping matrix) for next time step. 

        indexDrawdown = drawdownRun(:,t) >= drawdownLimit(n);
        if t < length(time)
            temp = repmat(QTheisRun(t,indexDrawdown)/1.2,length(time) - t, 1);
            QTheisRun(t+1:end,indexDrawdown) = temp;
        end
        

    end
    
    QRestricted(:,:,n) = QTheisRun;
    unmetDemand(:,:,n) = QInput - QRestricted(:,:,n); 
    drawdown(:,:,n) = drawdownRun;
    
   
end
    




