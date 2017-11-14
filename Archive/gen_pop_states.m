function [ s_pop, pop_M, T_pop ] = gen_pop_states( popParam)
% Generates population state space based on initial population,
% discretization, range of growth rates


% Generate possible growth rates
growthRates = popParam.min_growth: popParam.discrete_step_growth : popParam.max_growth;
growth_M = length(growthRates);

% Generate population states
s_pop = popParam.pop_initial: popParam.discrete_step_pop: popParam.pop_max;
pop_M = length(s_pop);

% Calculate next population state for valid growth states
growthMat = ones(1, growth_M) .* growthRates;
growthMat = repmat(growthMat,pop_M,1);
currentPop = ones(pop_M, 1) .* s_pop';
currentPop = repmat(currentPop, 1, growth_M);
nextPopUnrounded = currentPop .* ( ones(pop_M,growth_M) + growthMat );
nextPop = round2x(nextPopUnrounded,s_pop);

% Calculate transition probabilities
T_pop = zeros(pop_M);
for i = 1:pop_M
    index = ismember(s_pop, nextPop(i,:));
    countIndex = sum(index);
    prob = 1/countIndex;
    T_pop(i,index) = prob; 
end

