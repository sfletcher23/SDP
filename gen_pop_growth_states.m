function [s_pop, s_growth, T_growth, nextPop, max_end_pop] = gen_pop_growth_states(popParam, N)

% gen_pop_growth_states: Generates state space and transition informaiton 
% for input paramters about population growth

    % Outputs
    % s_pop: a vector with the possible population states
    % s_growth: a vector with the possible growth states
    % T_growth: a transition matrix for the growth rate
    % nextPop: a matrix [num pop states x num growth states] in which each
    % element is the population in t+1 for row and column indexes
    % reflecting population and growth rate in t
    
% Inputs
pop_initial = popParam.pop_initial;
min_growth = popParam.min_growth;
max_growth = popParam.max_growth;
max_growth_delta = popParam.max_growth_delta;
discrete_step_pop = popParam.discrete_step_pop;
discrete_step_growth = popParam.discrete_step_growth;
    
    
% Calculate upper and lower bound after N years
max_end_pop = pop_initial * (1 + max_growth) ^N;
min_end_pop = pop_initial * (1 + min_growth) ^N;
min_overall_pop = min(pop_initial, min_end_pop);
max_overall_pop = max(pop_initial, max_end_pop);

% Calculate discrete state space
s_pop = floor(min_overall_pop) : discrete_step_pop : ceil(max_overall_pop);
s_growth = min_growth : discrete_step_growth : max_growth;
growth_M = length(s_growth);
pop_M = length(s_pop);

% Calculate growth rate transition probabilities
% Assume uniform discrete distribution over valid range; can make
% triangular or similar in the future
T_growth = zeros(length(s_growth));
for i = 1:growth_M
    mid = s_growth(i);
    lb = mid - max_growth_delta;
    ub = mid + max_growth_delta;
    index = lb <= s_growth & s_growth <= ub; 
    prob_val = 1/sum(index);
    probArray = zeros(1,length(s_growth));
    probArray(index) = prob_val;
    T_growth(i,:) = probArray;
end

% Test Valid Prob Matrix
margin =  discrete_step_growth/1000;
err = sum(T_growth,2) - 1;
if sum(err > margin) 
    error('Invalid T_growth matrix')
end

% Calculate next population state for valid growth states
growthMat = ones(1, growth_M) .* s_growth;
growthMat = repmat(growthMat,pop_M,1);
currentPop = ones(pop_M, 1) .* s_pop';
currentPop = repmat(currentPop, 1, growth_M);
nextPopUnrounded = currentPop .* ( ones(pop_M,growth_M) + growthMat );
nextPop = round2x(nextPopUnrounded,s_pop);

% % Check discretization error
% err = abs(nextPop - nextPopUnrounded);
% margin = discrete_step_pop;
% if sum(err > margin)
%     error('Invaild Next Population State Matrix')
% end
% Throws error because when at high states and have high error, limited.
% Should never get there though. Put a test in somewhere else that the
% population never gets that big.

