%% Test parfor



% State definitions
s_expand = 1:2;
exp_M = length(s_expand); % Desalination expanded = 2

V = NaN(gw_M, exp_M, growth_M, N+1);
X1 = NaN(gw_M, exp_M, growth_M, N+1);
X2 = NaN(gw_M, exp_M, growth_M, N+1);

% Terminal period
X1(:,:,:,N+1) = zeros(gw_M, exp_M,  growth_M, 1);
X2(:,:,:,N+1) = zeros(gw_M, exp_M,  growth_M, 1);
V(:,:,:,N+1) = zeros(gw_M, exp_M,  growth_M, 1);



% Loop over all time periods
for t = linspace(N,1,N)
    
    % Calculate range of possible population states this time step
    s_pop_thisPeriod = pop_states_this_period(s_pop, t, nextPop);
    num_pop_thisPeriod = length(s_pop_thisPeriod);
        % Check that subset of total pop states
        isSubset = ismembertol(s_pop_thisPeriod,s_pop, 1E-4);
        test = sum(~isSubset);
        if test > 0
            error('Pop growth state set this period invalid')
        end
    
%     % Calculate nextV    
%     nextV = V(:,:,:,t+1);
    
    % Get T S index for this period
    index_T_S_samples_thisPierod = index_T_S_samples(:,t);
    
    t_now = t;
    
    X1_temp = X1(:,:,:,t);
    % Loop over all states
    % Loop over groundwater state: 1 is depleted, M1 is full
    parfor index_s1 = 1:gw_M
        s1 = s_gw(index_s1);
        
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for index_s2 = 1:exp_M
            s2 = s_expand(index_s2);
            
                
                % Loop over growth state
                for index_s4 = 1:growth_M
                    s4 = s_growth(index_s4);
                
                    bestV = Inf;
                    bestX = [0 0];  % Groundwater action and expansion action

                    % Update available actions based on whether gw depleted
                    if s1 == max(s_gw) 
                        a_gw = a_gw_depleted;
                    else
                        a_gw = a_gw_available;
                    end
 
                    % Update available actions based on whether expansion available
                    if s2 < 2
                        a_expand = a_expand_available;
                    else
                        a_expand = a_expand_unavailable;
                    end
                  
                    num_a_gw = length(a_gw);
                    num_a_expand = length(a_expand);
                    
                    % Loop over all actions
                    % Loop over groundwater pumping action
                    for index_a1 = 1:num_a_gw
                        a1 = a_gw(index_a1);
                        
                        % Loop over expansion action: 1 is do not expand, 2 is expand
                        for index_a2 = 1:num_a_expand
                            a2 = a_expand(index_a2); 
                        
                           bestX1 = a1;
                           bestX2 = a2;
                        end
                    end
                    
                    bestV = 3;
                    bestX1 = 1;
                    bestX2 = 2;
                    
                    % Save best value and action for current state
                    V(index_s1, index_s2,  index_s4, t_now) = bestV;
                    X1_temp = bestX1;
                    %X2(index_s1, index_s2,  index_s4, t_now) = bestX2;

                end
            end
        end
        X1(:,:,:,N+1) = X1_temp;
end



