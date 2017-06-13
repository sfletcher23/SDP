%% Current

T = cell(1,sizeA);
T_gw_exp = T;
for i = 1:sizeA
    T_gw_exp{i} = sparse(gw_M * exp_M, gw_M * exp_M);
    a = A(i,1);
    b = A(i,2); 
    range = gw_M*b+1:gw_M*(b+1);    % select which range is fille depending on exp state
    T_gw_exp{i}(gw_M+1:gw_M*2,gw_M+1:gw_M*2) = T_gw{a+1};  % 4th quadrant always filled
    T_gw_exp{i}(1:gw_M, range) = T_gw{a+1}; % Fill relevant 2nd range
end

clear T_gw T_exp

for i = 1:sizeA
    tempTCell = cell(pop_M);
    for j = 1:pop_M
        for k = 1:pop_M
            tempTCell{j,k} = T_gw_exp{i} * T_pop(j,k);
        end
    end
    T{i} = cell2mat(tempTCell);
end
%% Brute Force transition mat calculation

T = zeros(S);

    
% Calculate demand
% Population based on S index
indexPop = ceil((1:S) ./ (gw_M * exp_M));
pop = s_pop(indexPop,1);
dmd = water.demandPerCapita * 365/1000  ...    % m^3/p/year
    * pop * 1E6 ...   % p
    * water.demandFraction;
    



for i = 1:S

    
    
end
