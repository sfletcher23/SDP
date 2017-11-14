function [ s_pop, pop_M, T_pop ] = gen_pop_states_fromxls()

% Read pop and growth from table
s_pop = xlsread('population_states','PopStates');
pop_values = min(s_pop):.25:max(s_pop);
growth_values = 3:.5:6.5;
[pop_M,~] = size(s_pop);

% Calculate transition vector
T_pop = zeros(pop_M);
for i = 1:pop_M
    currentPop = s_pop(i,1);
    currentGrowth = s_pop(i,2);
    nextPopUnrounded = currentPop * (1 + currentGrowth/100);
    nextPop = round2x(nextPopUnrounded, pop_values);
    indexNextPop = find(nextPop  == s_pop(:,1));
    indexNextPopCurrentGrowth = find(s_pop(indexNextPop,2) == currentGrowth) + indexNextPop(1) -1;
    % Check current growth is available
    if isempty(indexNextPopCurrentGrowth)
        indexNextPopCurrentGrowth = max(indexNextPop);
    end
    T_pop(i,indexNextPopCurrentGrowth) = 1/2;
    % Check if lower growth is available
    checkLow = find(indexNextPopCurrentGrowth - 1  == indexNextPop,1);
    if isempty(checkLow)
        T_pop(i,indexNextPopCurrentGrowth) = T_pop(i,indexNextPopCurrentGrowth) + 1/4;
    else
        T_pop(i,indexNextPopCurrentGrowth - 1) = 1/4;
    end
    % Check if upper growth is available
    checkHigh = find(indexNextPopCurrentGrowth + 1  == indexNextPop,1);
    if isempty(checkHigh)
        T_pop(i,indexNextPopCurrentGrowth) = T_pop(i,indexNextPopCurrentGrowth) + 1/4;
    else
        T_pop(i,indexNextPopCurrentGrowth + 1) = 1/4;
    end
end

%Check T
rowSums = sum(T_pop,2);
indexCheck = ~(rowSums == 1);
if sum(indexCheck) > 0 
    error('Invalid T_pop')
end



end