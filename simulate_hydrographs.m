%% Simulate hydrograph
if false
load('sample_nn_inputstargets.mat')
end

%% Get nn function

nnNumber = 17182;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
adjustOutput = true;
wellIndex = 108;
saveOn = true;

%% Get modflow, nn, and sdp hydrograph estimates and plot: One well one run
if false

% Sample run number and well number to plot
numRuns = length(hk);
i = randsample(numRuns,1);


% Get row index range of time series corresponding to sample
indexMin = (i-1)*365*30 + 1;
indexMax = i*365*30;

% Get estimates from nn and modflow output
xsample = x(:,indexMin:indexMax);
y_nn = netscript(xsample, adjustOutput);
y_modflow = t(:, indexMin:indexMax);
y_nn = y_nn(wellIndex,:);
y_modflow = y_modflow(wellIndex,:);

% Get estimates using SDP
s_gw = 0:200;
N = 30;
gw_state = zeros(1,N);
% [K_samples_thisPeriod, S_samples_thisPeriod] = gen_param_dist('full_range', gwParam, 1, N);
K_samples_thisPeriod = hk(i);
S_samples_thisPeriod = sy(i);
for time = 1:N
    % Get transmat vector to next GW state 
    if time == 1
        gw_state_previous = 0;
    else
        gw_state_previous = gw_state(time-1);
    end
    T_current_gw = gw_transrow_nn(nnNumber, wellIndex, time, K_samples_thisPeriod, S_samples_thisPeriod, gw_state_previous, s_gw, adjustOutput);
    p = rand();
    index = find(p < cumsum(T_current_gw),1);
    gw_state(time) = s_gw(index);
end
y_sdp = 200 - gw_state;

% Plot
figure;
plot([1:30*365]/365, y_modflow)
hold on
plot([1:30*365]/365, y_nn)
plot(1:30, y_sdp)
legend('modflow', 'nn', 'sdp')

% clear y_nn y_modflow y_sdp gw_state gw_state_previous T_current_gw K_samples_thisPeriod S_samples_thisPeriod time index

end
%% Simulate sdp hydrographs from multiple using data updating

sampleSize = 10000;
runs = 20;
N = 30;
s_gw = 0:200;
gw_state = zeros(runs,N);
numSampUsed = zeros(runs,N);
[K_samples, S_samples] = gen_param_dist('full_range', sampleSize, 1, N);
parfor i = 1:runs
    disp(num2str(i))
    tempGwState = zeros([1 N]);
    for t = 1:N
        % Get transmat vector to next GW state 
        gw_state_current = tempGwState(t);
        [T_current_gw, numSampUsedNow] = gw_transrow_nn(nnNumber, wellIndex, t, K_samples, S_samples, gw_state_current, s_gw, adjustOutput);
        numSampUsed(i,t) = numSampUsedNow;
        p = rand();
        index = find(p < cumsum(T_current_gw),1);
        tempGwState(t) = s_gw(index);
    end
    gw_state(i,:) = tempGwState;
end
y_sdp = 200 - gw_state;

figure;
plot(1:30, y_sdp)
ylim([0 200])

if saveOn
    datetime=datestr(now);  
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save(strcat('hydrograph_sim_', datetime), 'y_sdp')
end
