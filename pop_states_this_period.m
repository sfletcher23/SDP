function [s_pop_thisPeriod] = pop_states_this_period(s_pop, t, nextPop)

% Finds the range of pop values available during this time period. It can
% be higher than the growth rate would suggest because of the
% discretization. Therefore, we cycle through nextPop at the maximum
% transision t times to find the max.

startingPop = min(s_pop);
for i = 1:t+1
    startingPopIndex = find(s_pop == startingPop);
    maxPopNextPeriod = nextPop(startingPopIndex,end);
    startingPop = maxPopNextPeriod;
end

maxPopIndex = find(s_pop == maxPopNextPeriod);
s_pop_thisPeriod = s_pop(1:maxPopIndex);

end