function[T_gw] = gw_transrow_nn(nnNumber,wellIndex, t, K_samples_thisPeriod, S_samples_thisPeriod, s1, s_gw ) 

% Calculates drawdown between time t-1 and time t predicted by the neural
% net indicated for each of the T and S samples indicated. 
% Uses those drawdown estimates to develop probability distribution for
% next groundwater state

netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 

% Generate inputs to neural net function for t and t-1
% Input order: hk, sy, time  Input size: 3 x numSamples*2

timeStepSize = 365;

[~, numSamples] = size(K_samples_thisPeriod);
time = repmat(t*timeStepSize, [1 numSamples]);
x = [K_samples_thisPeriod; S_samples_thisPeriod; time];
head_t = netscript(x);
head_t = head_t(wellIndex,:);

if t>1
    time = repmat(t-1, [1 numSamples]);
    x = [K_samples_thisPeriod; S_samples_thisPeriod; time];
    head_t_previous = netscript(x);
    head_t_previous = head_t_previous(wellIndex,:);
else
    head_t_previous = repmat(200, [1 numSamples]);
end
    
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
