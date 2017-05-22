function [T_S_pair_cdf] = gen_param_dist(T_S_pairs, t)
% Generate joint CDF across the aquifer parameters
% For now, assume uniform distribution with 10% cut off the top of the
% range of both T and S each time step

T = T_S_pairs(:,1);
S = T_S_pairs(:,2);
T_upper = max(T);
T_lower = min(T);
T_median = median(T);
S_upper = max(S);
S_lower = min(S);
S_median = median(S);
[numPairs, ~] = size(T_S_pairs);
T_S_pair_pmf = zeros(numPairs,1);

% Cutoff top 10% of T and S each time, assume uniform.
T_90 = T_upper * .95^t;
S_90 = S_upper * .95^t;
indexNonZeroT = T_S_pairs(:,1) < T_90;
indexNonZeroS = T_S_pairs(:,2) < S_90;
indexNonZero = indexNonZeroT .* indexNonZeroS;
countNonZero = sum(indexNonZero);
uniProb = 1/countNonZero;
T_S_pair_pmf(logical(indexNonZero)) = uniProb; 
T_S_pair_cdf = cumsum(T_S_pair_pmf);


end