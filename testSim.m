%% Test simulation results

% GW state should monotonically decrease

for t = 2:length(state_gw)
    diff = state_gw(t) - state_gw(t-1);
    if diff > 0 
        error(['GW state not monotonically decreasing in year ' num2str(t)])
    end
end
