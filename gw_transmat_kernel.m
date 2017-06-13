function [T_gw] = gw_transmat_kernel(gw_supply, kernel, index_T_S_samples, s_gw, gw_M )

% Calculates the range of drawdown (dd_values) and the probability of
% those values (dd_prob) for a given distribution of possible parameters
% and a pumping volume

% For now, using only 1st observation well and pumping well

% kernel is a [numObserve x numTime x numParameterValues x numPumpWells] matrix


% Calculate drawdown for each paramter sample using kernel functions
kernel_sampled = kernel(:,2,index_T_S_samples,:);
drawdown = kernel2drawdown(gw_supply, kernel_sampled, 2 );
drawdown = reshape(drawdown, [1, numel(drawdown)]);

% Calculate next state
next_gw_state = repmat(s_gw,[1 length(drawdown)]) + repmat(drawdown, [gw_M, 1]);
rounded_next_gw_state = round2x(next_gw_state, s_gw);

% Calculate transition mat
T_gw = zeros(gw_M, gw_M);
for i = 1:gw_M
T_gw(i,:) = histcounts(rounded_next_gw_state(i,:),  [s_gw' s_gw(end)+1], 'Normalization', 'probability');
end

% Test valid prob distribution
margin = 1E-4;
err = abs(sum(T_gw,2) - 1);
if err > margin
    error('invalid probability distribution for T_gw')
end



end

