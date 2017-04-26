function [T_gw] = gw_transrow_kernel(gw_supply, kernel, index_T_S_samples, t, s1, s_gw )

% Calculates the range of drawdown (dd_values) and the probability of
% those values (dd_prob) for a given distribution of possible parameters
% and a pumping volume

% For now, using only 1st observation well and pumping well

% kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix

% Calculate drawdown for each paramter sample using kernel functions
Q = gw_supply;
kernel_sampled = kernel(:,t,index_T_S_samples,:);
drawdown = kernel2drawdown(Q, kernel_sampled, t );

% Calculate next state
next_s1 = s1 + drawdown;
rounded_next_s1 = round2x(next_s1, s_gw);

% Calculate transition probability row
T_gw = histcounts(rounded_next_s1,  [s_gw s_gw(end)+1], 'Normalization', 'probability');

% Test valid prob distribution
margin = 1E-4;
err = abs(sum(T_gw) - 1);
if err > margin
    error('invalid probability distribution for T_gw')
end



end

