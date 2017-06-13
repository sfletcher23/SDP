%% Compute reward matrix
tic
R = zeros(S, sizeA);
for s1 = s_gw'
    index_s1 = find(s1 == s_gw);
    for s2 = s_expand'
        index_s2 = find(s2 == s_expand);
        for index_s3 = 1:pop_M
            pop = s_pop(index_s3,1);
            for a1 = gw_A'
                for a2 = exp_A'

                    % If this action is not available at this state, cost = Inf
                    if (s2 == 2 && a2 == 1) || (s1 == s_gw(end) && a1 == 1)
                        RV = 1e30;
                    else % Calculate cost and shortages
                        demand = get_demand(water, pop, water.demandFraction); % Calculate demand
                        [ RV, ~, ~, ~, ~, ~, ~] = getCost(a1, a2, s1, s2, water, demand, s_gw, gwParam, costParam);  % Can vectorize this fucntion later
                    end

                    % Calculate action index
                    indexAction = vectorIndex([a1+1 a2+1], A_vectors);   % +1 since indexing starts at 1, actions start at 0, maybe change action definitions?

                    % Calculate state vectors
                    indexState = vectorIndex([index_s1 index_s2 index_s3], {s_gw, s_expand, s_pop});

                    % Put cost in reward matrix
                    R(indexState, indexAction) = -RV;

                end
            end
        end
    end
end
toc

tic

% Vectorized version
RV = zeros(S, sizeA);

% Calculate demand
% Population based on S index
indexPop = ceil((1:S) ./ (gw_M * exp_M));
pop = s_pop(indexPop,1);
demand = water.demandPerCapita * 365/1000  ...    % m^3/p/year
    * pop * 1E6 ...   % p
    * water.demandFraction;

% Calculate supply
% Initial desal supply
supply = ones(size(RV)) * water.desal_capacity_initial;
% Add expansion supply
ind1 = mod(1:S,gw_M*2) > gw_M;
ind2 = mod(1:S,gw_M*2) == 0;
indexExpanded = ind1 | ind2;
supply(indexExpanded,:) = supply(indexExpanded,:) + water.desal_capacity_expansion; 
% Add pumping supply
gw_supply = zeros(size(RV));
gw_supply(:,1) = 0;  % no pump, no exp
gw_supply(:,2) = gwParam.pumpingRate;  % pump, no exp
gw_supply(:,3) = 0;  % no pump, exp
gw_supply(:,4) = gwParam.pumpingRate;  % pump, exp
supply = supply + gw_supply;

% Calculate shortage
shortage = max(0, demand - supply);

% Costs include shortage costs, expansion costs, and pumping costs
shortageCost = shortage * costParam.shortage_cost;
expansionCost = zeros(size(RV));
expansionCost(:,3:4) = costParam.expansion_cost;
pumpingCost = costParam.pumping_cost * gw_supply;
RV = -(shortageCost + expansionCost + pumpingCost);

% Infeasible actions
infeasibleCost = 1e30;
% Pumping with max drawdown
indexMaxDrawdown = (mod(1:S,gw_M) == 0);
RV(indexMaxDrawdown, 2) = -infeasibleCost;
RV(indexMaxDrawdown, 4) = -infeasibleCost;
% Expanding when already expanded
RV(indexExpanded,3:4) = -infeasibleCost;

toc
% Check R and RV
diff = R - RV;
err = 0.01;
index = diff > err;
sum(sum(index))


