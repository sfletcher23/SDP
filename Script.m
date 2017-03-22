%% Simple SDP GW Model

%% Definitions / Initialize

% Time period
N = 6;

% Number states
M1 = 5; % Groundawter
M2 = 2; % Desalination expanded = 2

% Actions: groundwater pumping high or low or none, expand desal
% a1 pumping actions: 0 no pumping, 1 low pumping, 2 high pumping
a1_available = [0 1 2];
a1_depleted = [0];
% a2 desal actions: 0 no expand, 1 expand
a2_available = [0 1];
a2_unavailable = [0];

% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = zeros(M1, M2, N+1);
X = zeros(M1, M2, N+1);

% Terminal period
X(:,:,N+1) = zeros(M1, M2, 1);
V(:,:,N+1) = zeros(M1, M2, 1);

%% Cost parameters

demand = 10;
desal_capacity = 5;
shortage_cost = 1;
expansion_cost = 4;

%% Transition Matrix

% T2: high pumping
T2 = zeros(M1,M1);
T2(1,1) = 1;
for i = 2:M1
    T2(i,i) = 0.5;
    T2(i,i-1) = 0.5;
end

% T1: low pumping
T1 = zeros(M1,M1);
T1(1,1) = 1;
for i = 2:M1
    T1(i,i) = 0.8;
    T1(i,i-1) = 0.2;
end

% T0: no pumping
T0 = zeros(M1,M1);
T0(1,1) = 1;
for i = 2:M1
    T0(i,i) = 1;
    T0(i,i-1) = 0;
end


%% Backwards Recursion

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Loop over all states
    % Loop over groundwater state: 1 is depleted, M1 is full
    for s1 = 1:M1       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for s2 = 1:M2

            bestV = 999999;
            bestX = [0 0];  % Groundwater action and expansion action

            % Update available actions based on whether GW depleted
            if s1 > 1 
                a1 = a1_available;
            else
                a1 = a1_depleted;
            end
            
            % Update available actions based on whether expansion available
            if s2 < 2
                a2 = a2_available;
            else
                a2 = a2_unavailable;
            end

            % Loop over all actions
            % Loop over groundwater pumping action
            for a1 = 1:length(a1)
                % Loop over expansion action: 1 is do not expand, 2 is expand
                for a2 = 1:length(a2)
                    
                    % Calculate cost this period
                    capacity = s1 - 1;   % capacity from GW 0 to 4
                    capacity = capacity + desal_capacity * s2; % desal if available
                    shortage = max(0, demand - capacity);
                    % Cost = shortage cost + expasion cost if expanded
                    cost = shortage * shortage_cost + expansion_cost * a2;
                    
                    
                    % Pick transmat row for GW based on action, current
                    % GW state
                    if a1 == 0 
                        Ta1 = T0(s1,:);
                    elseif a1 == 1
                        Ta1 = T1(s1,:);
                    else
                        Ta1 = T2(s1,:);
                    end
                    
                    % Choose V based on which desal state
                    if a2 == 2 || s2 == 2   % desal already expanded or will expand
                        nextV = V(:,2,t+1);
                    else
                        nextV = V(:,1,t+1);
                    end
                    
                   % Calculate expected future cost
                   expV = Ta1 * nextV;
                                            
                   % Check if best decision
                   checkV = cost + expV;
                   if checkV < bestV
                       bestV = checkV;
                       bestX = [a1 a2];
                   end
                   
                end
            end
            
            % Save best value and action for current state
            V(s1, s2, t) = bestV;
            X1(s1, s2, t) = bestX(1);
            X2(s1, s2, t) = bestX(2);
            
        end
    end
end

% Pumping doesn't actually contribute to reward here! So just picks first
% one every time.

