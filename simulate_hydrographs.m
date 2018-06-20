%% Simulate hydrograph
if true
load('sample_modflow_data.mat')
end

%% Get nn function

nnNumber = 54212;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
adjustOutput = true;
nstp = 100;
saveOn = true;
% GW Parameters
gwParam = struct;
gwParam.initialDrawdown = 0;
gwParam.sampleSize = 100;
gwParam.depthLimit = 300;
gwParam.pumpingRate = 640000 * 365;  % m^3/y
gwParam.otherPumpingRate = (970000 + 100000 - 640000) * 365;  % m^3/y    % From ADA water balance report 2016 estimates
gwParam.nnNumber = 54212;
gwParam.wellIndex = 108; % 68 is RR1, 108 is Shemess, 93 is royal garage
gwParam.exaggeratePumpCost = false;
gwParam.enforceLimit = false;
gwParam.pumpingSubsidy = true;
gwParam.infoScenario = 'full_range';
gwParam.TgwLoadName = 'T_gw';
gwParam.likelihoodfct = 'normal';
gwParam.llhstddev = 10;
gwParam.startingHead = 337.143;
gwParam.nstp = 100;
%% Get modflow, nn, and sdp hydrograph estimates and plot: One well one run
if true
    
% Sample run number and well number to plot
numRuns = length(hk);
i = randsample(numRuns,1);

N = 30;

% Get row index range of time series corresponding to sample
indexMin = (i-1)*nstp*N + 1;
indexMax = i*nstp*N;

% Get estimates from nn and modflow output
xsample = x(:,indexMin:indexMax);
y_nn = netscript(xsample, gwParam);
y_modflow = t(:, indexMin:indexMax);
time_nn = xsample(3,:);
time_sdp = 1:N;
indexTime = ismember(round(time_nn/365,2),time_sdp);
time_nn = time_nn(indexTime)/365;
y_nn = y_nn(indexTime);
y_modflow = y_modflow(indexTime);


% Get estimates using SDP
[s_gw, gw_M] = gen_water_growth_states(gwParam);
gw_state = zeros(1,N);
% [K_samples_thisPeriod, S_samples_thisPeriod] = gen_param_dist('full_range', gwParam, 1, N);
K_samples_thisPeriod = exp(hk(i));
S_samples_thisPeriod = exp(ss(i));
% K_samples_thisPeriod = 0.6701;
% S_samples_thisPeriod = 0.2062;

for time = 1:N-1
    % Get transmat vector to next GW state 
    if time == 1
        gw_state_current = 0;
    else
        gw_state_current = gw_state(time);
    end
    [T_current_gw, numRelevantSamples, stateInfeasible, indexAbove, indexBelow, indexRelevantSamples, drawdown] = gw_transrow_nn(gwParam, time, K_samples_thisPeriod, S_samples_thisPeriod, gw_state_current, s_gw);
    p = rand();
    index = find(p < cumsum(T_current_gw),1);
    gw_state(time+1) = s_gw(index);
end
y_sdp = gwParam.startingHead - gw_state;
%y_sdp = y_sdp(2:end);

% Plot
figure;
plot(time_nn + 1, y_modflow)
hold on
plot(time_nn + 1, gwParam.startingHead - y_nn)
time_sdp = 1:N;
plot(time_sdp, y_sdp)
legend('modflow', 'nn', 'sdp')
ylim([-650 400])
% clear y_nn y_modflow y_sdp gw_state gw_state_previous T_current_gw K_samples_thisPeriod S_samples_thisPeriod time index

end

%% Error histogram

if true
    runs = round(1000/3);
    N = 30;
    error_nn = zeros([length(runs) N-1]);
    error_sdp = zeros([length(runs) N]);
    error_sdp_nn = zeros([length(runs) N]);
    y_nn = zeros([length(runs), N]);
    y_modflow = zeros([length(runs), N]);
    
    for j = 1:runs
        
        i = randsample(numRuns,1);

        % Get row index range of time series corresponding to sample
        indexMin = (i-1)*nstp*N + 1;
        indexMax = i*nstp*N;

        % Get estimates from nn and modflow output
        xsample = x(:,indexMin:indexMax);
        y_nn_temp = netscript(xsample, gwParam);
        y_modflow_temp = t(:, indexMin:indexMax);
        time_nn = xsample(3,:);
        time_sdp = 0:N;
        indexTime = ismember(round(time_nn/365,2),time_sdp); 
        y_nn(j,:) =  y_nn_temp(indexTime);
        y_modflow(j,:) = y_modflow_temp(indexTime);
        

%         % Get estimates using SDP
%         [s_gw, gw_M] = gen_water_growth_states(gwParam);
%         gw_state = zeros(1,N+1);
%         % [K_samples_thisPeriod, S_samples_thisPeriod] = gen_param_dist('full_range', gwParam, 1, N);
%         K_samples_thisPeriod = hk(i);
%         S_samples_thisPeriod = ss(i);
%         % K_samples_thisPeriod = 0.6701;
%         % S_samples_thisPeriod = 0.2062;
%         for time = 1:N
%             % Get transmat vector to next GW state 
%             if time == 1
%                 gw_state_current = 0;
%             else
%                 gw_state_current = gw_state(time);
%             end
%             [T_current_gw, numRelevantSamples, stateInfeasible, indexAbove, indexBelow, indexRelevantSamples, drawdown] = gw_transrow_nn(gwParam, time, K_samples_thisPeriod, S_samples_thisPeriod, gw_state_current, s_gw, adjustOutput);
%             p = rand();
%             index = find(p < cumsum(T_current_gw),1);
%             gw_state(time+1) = s_gw(index);
%         end
%         y_sdp = gwParam.startingHead - gw_state;
%         y_sdp = y_sdp(2:end);
%         y_nn = gwParam.startingHead - y_nn;
        
        
%         error_sdp(j,:) = y_sdp - y_modflow;
%         error_sdp_nn(j,:) = y_sdp - y_nn;

    end
    
    dd_nn = diff(gwParam.startingHead - y_nn,1, 2);
    dd_modflow = diff(y_modflow,1, 2);
    error_nn = dd_nn - dd_modflow;
    index_above_depth = y_modflow > (gwParam.startingHead - gwParam.depthLimit);
    error_above_depth = error_nn(index_above_depth(:,2:end));
    
    figure;
%     subplot(1,3,1) 
    hist(error_above_depth,50,'k');
    h = gca;
    h.ColorOrder = repmat([1 1 1], [7, 1]);
    set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
    xlabel('ANN drawdown - MODFLOW drawdown [m]')
    ylabel('instances')
    rmse = sqrt(mean(reshape(error_above_depth, 1, []) .^ 2));
    title(strcat('Drawdown Error from ANN: RMSE= ', num2str(rmse)))
    xlim([-20 20])
    
%     subplot(1,3,2)
%     hist(error_sdp_nn,10,'k');
%     h = gca;
%     h.ColorOrder = repmat([1 1 1], [7, 1]);
%     set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
%     xlabel('SDP head - ANN head')
%     ylabel('instances')
%     rmse = sqrt(mean(reshape(error_sdp_nn, 1, []) .^ 2));
%     title(strcat('RMSE: ', num2str(rmse)))
%     xlim([-20 20])
%     
%     subplot(1,3,3)
%     hist(error_sdp,10,'k');
%     h = gca;
%     h.ColorOrder = repmat([1 1 1], [7, 1]);
%     set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
%     xlabel('SDP head - MODFLOW head')
%     ylabel('instances')
%     rmse = sqrt(mean(reshape(error_sdp,1, []) .^ 2));
%     title(strcat('RMSE: ', num2str(rmse)))
%     xlim([-20 20])
    
    
end

%% Simulate sdp hydrographs from multiple using data updating
if false

sampleSize = 1000;
runs = 20;
N = 30;
[s_gw, gw_M] = gen_water_growth_states(gwParam);
gw_state = zeros(runs,N+1);
numSampUsed = zeros(runs,N);
drawdown = cell(runs,N);
[K_samples, S_samples] = gen_param_dist(sampleSize);

figure;
for i = 1:runs
    disp(num2str(i))
    tempGwState = zeros([1 N+1]);
    tempNumSamples = zeros([1 N]);
    tempDrawdown = zeros([1 N]);
    for time = 1:N
        % Get transmat vector to next GW state 
        gw_state_current = tempGwState(time);
        index_s1 = find(gw_state_current == s_gw);
        [T_current_gw, numSampUsedNow, ~, ~, ~, ~, tempDrawdown] = gw_transrow_nn(gwParam, time, K_samples, S_samples, gw_state_current, s_gw);
        %T_current_gw = T_gw_all(:,index_s1,time);
        tempNumSamples(time) = numSampUsedNow;
        p = rand();
        index = find(p < cumsum(T_current_gw),1);
        tempGwState(time+1) = s_gw(index);
    end
    gw_state(i,:) = tempGwState;
    %numSampUsed(i,:) = tempNumSamples;
end
y_sdp = gwParam.startingHead - gw_state;

hold on
plot(0:30, y_sdp)
ylim([-650 350])

end

if false
    datetime=datestr(now);  
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save(strcat('hydrograph_sim_', datetime), 'y_sdp')
end
