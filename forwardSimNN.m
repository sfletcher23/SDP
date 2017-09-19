

%% Gen samples
nnNumber = 17182;
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 
adjustOutput = true;
wellNum = 108;
infoScenario = 'full_range';
sampleSize = 10000;
t = 1;
N = 30;
[K_samples, S_samples] = gen_param_dist(infoScenario, sampleSize, t, N);
%% transrow
s_gw = 0:200;
for t = 1:N
    for s1 = s_gw
    [T_gw, numRelevantSamples, stateInfeasible, indexAbove, indexBelow, indexRelevant] = gw_transrow_nn(nnNumber,wellIndex, t, K_samples, S_samples, s1, s_gw, adjustOutput); 
    end
end

%% Forward simulate NN to try to get mid range drawdown for a well

time = 1:365*30;
n = length(time);
K_lower = 0.4; % [ m^2/day]
K_upper = 2.5;
K_mean = 1.170;
K_sigma = 0.56;
S_lower = 0.02;
S_upper = 3e-1;
S_mean = 0.13;
figure;
for i = 1:100
    hk = K_samples(i);
    sy = S_samples(i);
    x = [ones(1,n)*hk; ones(1,n)*sy; time];
    y = netscript(x, true);
    y = y(wellNum,:);
    hold on
    plot(time/365, y)
end
line([30 30], [0 200])
ylim([0 200])

%% 

