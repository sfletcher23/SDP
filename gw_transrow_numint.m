function[T_gw] = gw_transrow_numint(gwParam, s1, s_gw, drawdown ) 

% If reach depth limit, transition to -1
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

% Calculate next state
next_s1 = s1 + drawdown;
rounded_next_s1 = round2x(next_s1, s_gw);

% Calculate transition probability row
T_gw = histcounts(rounded_next_s1,  [s_gw(1:end):s_gw(end)] + 0.1, 'Normalization', 'probability');
T_gw = [0 T_gw];

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
    save(strcat('t_gw_error_', datetime), 'nnNumber' , 't', 'K_samples', 'S_samples', 's1', 's_gw');
    error(strcat('Invalid probability distribution for T_gw'))
end

end
