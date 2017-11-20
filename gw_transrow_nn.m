function[T_gw, numRelevantSamples, stateInfeasible, indexAbove, indexBelow, sampleProb, drawdown] = ...
    gw_transrow_nn(gwParam, t, K_samples, S_samples, s1, s_gw, adjustOutput ) 

% Calculates drawdown between time t-1 and time t predicted by the neural
% net indicated for each of the T and S samples indicated. 
% Uses those drawdown estimates to develop probability distribution for
% next groundwater state

nnNumber = gwParam.nnNumber;
nstp = gwParam.nstp;

stateInfeasible = false;
indexAbove = [];
indexBelow = [];
sampleProb = [];
drawdown = [];

K_samples = log(K_samples);
S_samples = log(S_samples);

% If at max drawdown, stay at max drawdown
% if s1 == 200
%     T_gw = zeros(1,length(s_gw));
%     T_gw(end) = 1;
%     numRelevantSamples = -888;
%     return
% end

% If there is a nonzero depth limit, transition to -1
if gwParam.depthLimit
    if s1 >= gwParam.depthLimit
        T_gw = zeros(1,length(s_gw));
        T_gw(1) = 1;
        numRelevantSamples = -7;
        return
    end
end

% If stopped pumping, stay at stopped pumping
if s1 == -1
    T_gw = zeros(1,length(s_gw));
    T_gw(1) = 1;
    numRelevantSamples = -99;
    return
end
 
% Get neural net script
netname = strcat('myNeuralNetworkFunction_', num2str(nnNumber));
netscript = str2func(netname); 

% Generate inputs to neural net function for t and t+1
% Input order: hk, sy, time  Input size: 3 x numSamples*2

[~, numSamples] = size(K_samples);

drawdown_t_current = zeros(1, length(K_samples));
drawdown_t_next = zeros(1, length(K_samples));
for i = 1:length(K_samples)
    x = [repmat(K_samples(i),[1,t+1]); repmat(S_samples(i),[1,t+1]); 0:365:365*(t) ];
    tempHead = netscript(x, gwParam);
%     if t == 1
%         drawdown_t_current(i) = 0;
%     else
%         drawdown_t_current(i) = tempHead(end-1);
%     end
    drawdown_t_current(i) = tempHead(end-1);
    drawdown_t_next(i) = tempHead(end);
end

if strcmp(gwParam.likelihoodfct, 'uniform') % Need to fix this with new drawdon approach if use again
    margin = gwParam.llhstddev; 
    indexRelevantSamples = abs(drawdown_t_current - s1) < margin;
    numRelevantSamples = sum(indexRelevantSamples);

    if numRelevantSamples == 0
        warning(strcat('infeasible groundwater state : t=',num2str(t), ', drawdown = ', num2str(s1), ...
            ', min nn head =', num2str(min(drawdown_t_current)), ', max nn head =', num2str(max(drawdown_t_current)) ));
        stateInfeasible = true;
        % Use closest sample even if outside error marign
        [~, bestIndex] = min(abs(drawdown_t_current - s1));
        indexRelevantSamples = false([1 numSamples]);
        indexRelevantSamples(bestIndex) = 1;
        indexRelevantSamples = logical(indexRelevantSamples);
        numRelevantSamples = 1;
        % Find above and below samples for analysis
        [~, indexAbove] = min(drawdown_t_current -s1);
        [~, indexBelow] = min(s1 - drawdown_t_current);
    end


    % In first period, accept all vlaues
%     if t == 1 
%         indexRelevantSamples = true([1 numSamples]);
%     end

    % Update samples from previous period to include only relevant samples
    drawdown_t_current = drawdown_t_current(indexRelevantSamples);
    drawdown_t_next = drawdown_t_next(indexRelevantSamples);

    % Calculate drawdown between t+1 and t
    drawdown =  drawdown_t_current - drawdown_t_next;
    indexNeg = drawdown < 0;
    drawdown(indexNeg) = 0;

    % Calculate next state
    next_s1 = s1 + drawdown;
    rounded_next_s1 = round2x(next_s1, s_gw);

    % Calculate transition probability row
    T_gw = histcounts(rounded_next_s1,  [s_gw(1:end):s_gw(end)] + 0.1, 'Normalization', 'probability');
    T_gw = [0 T_gw];


elseif strcmp(gwParam.likelihoodfct, 'normal')

    % Calculate probability for each parameter sample    
    u = drawdown_t_current;
    x = s1; 
    sampleProb = normpdf(x, u, gwParam.llhstddev);
    numRelevantSamples = sum(sampleProb > 0.001);
    
    % Calculate next state for each parameter sample
    drawdown =  drawdown_t_next - drawdown_t_current;
    indexNeg = drawdown < 0;
    drawdown(indexNeg) = 0;
    next_s1 = s1 + drawdown;
    rounded_next_s1 = round2x(next_s1, s_gw);
    
    % Sum probabilites for each possible next state
    zeropstates = s_gw(~ismember(s_gw,rounded_next_s1));
    [a,~,c] = unique([rounded_next_s1 zeropstates]);
    prob = [sampleProb zeros(1,length(zeropstates))];     
    out = [a', accumarray(c,prob)];
    sumout = sum(out(:,2));
    T_gw = out(:,2)'/sumout;
    sampleProb = sampleProb/sumout;
    
end

% Find states above depth limit with positive prob and switch to absorbing state
if gwParam.depthLimit
    indexAboveLimit = find(s_gw >= gwParam.depthLimit & T_gw > 0);
    sumAboveLimit = sum(T_gw(indexAboveLimit));
    T_gw(1) = sumAboveLimit;
    T_gw(indexAboveLimit) = 0;
end

% Test valid prob distribution
margin = 1E-4;
err = abs(sum(T_gw) - 1);
if err > margin
    datetime=datestr(now);  
    datetime=strrep(datetime,':','_'); %Replace colon with underscore
    datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
    datetime=strrep(datetime,' ','_');%Replace space with underscore
    save(strcat('t_gw_error_', datetime), 'nnNumber' , 't', 'K_samples', 'S_samples', 's1', 's_gw', 'adjustOutput');
    error(strcat('Invalid probability distribution for T_gw'))
end
