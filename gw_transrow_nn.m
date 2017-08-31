function[T_gw, numRelevantSamples, stateInfeasible, indexAbove, indexBelow] = gw_transrow_nn(nnNumber,wellIndex, t, K_samples_thisPeriod, S_samples_thisPeriod, s1, s_gw, adjustOutput ) 

% Calculates drawdown between time t-1 and time t predicted by the neural
% net indicated for each of the T and S samples indicated. 
% Uses those drawdown estimates to develop probability distribution for
% next groundwater state

stateInfeasible = false;
indexAbove = [];
indexBelow = [];

% If at max drawdown, stay at max drawdown
if s1 == s_gw(end)
    T_gw = zeros(1,length(s_gw));
    T_gw(end) = 1;
    numRelevantSamples = -999;
    return
end

% Get neural net script
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 

% Generate inputs to neural net function for t and t-1
% Input order: hk, sy, time  Input size: 3 x numSamples*2

timeStepSize = 365;
[~, numSamples] = size(K_samples_thisPeriod);
time = repmat(t*timeStepSize, [1 numSamples]);

% Find only samples close to current state: compare s1 and head at t-1
if t>1
    time = repmat((t-1)*timeStepSize, [1 numSamples]);
    x = [K_samples_thisPeriod; S_samples_thisPeriod; time];
    head_t_previous = netscript(x, adjustOutput);
    head_t_previous = head_t_previous(wellIndex,:);
else 
    head_t_previous = repmat(200, [1 numSamples]);
end
margin = 5; 
indexRelevantSamples = abs(head_t_previous - (200 -s1)) < margin;
numRelevantSamples = sum(indexRelevantSamples);
if numRelevantSamples == 0
    warning(strcat('infeasible groundwater state : t=',num2str(t), ', 200-s1 = ', num2str(200-s1), ...
        ', min previous head =', num2str(min(head_t_previous)), ', max previous head =', num2str(max(head_t_previous)) ));
    stateInfeasible = true;
    % Use closest sample even if outside error marign
    [~, bestIndex] = min(abs(head_t_previous - (200 -s1)));
    indexRelevantSamples = zeros([1 numSamples]);
    indexRelevantSamples(bestIndex) = 1;
    indexRelevantSamples = logical(indexRelevantSamples);
    numRelevantSamples = 1;
    % Find above and below samples for analysis
    [~, indexAbove] = min(head_t_previous - (200 -s1));
    [~, indexBelow] = min((200 -s1) - head_t_previous);
end

% Update samples from previous period to include only relevant samples
head_t_previous = head_t_previous(indexRelevantSamples);

% Get head estimates for next period from relevant samples
time = repmat(t*timeStepSize, [1 numRelevantSamples]);
x = [K_samples_thisPeriod(indexRelevantSamples); S_samples_thisPeriod(indexRelevantSamples); time];
head_t = netscript(x, adjustOutput);
head_t = head_t(wellIndex,:);
    
% Calculate drawdown between t and t-1
drawdown = head_t_previous - head_t;

% Calculate next state
next_s1 = s1 + drawdown;
rounded_next_s1 = round2x(next_s1, s_gw);

% Calculate transition probability row
T_gw = histcounts(rounded_next_s1,  [s_gw s_gw(end)+1], 'Normalization', 'probability');

% Test valid prob distribution
margin = 1E-4;
err = abs(sum(T_gw) - 1);
if err > margin
    error('invalid probability distribution for T_gw')
end
