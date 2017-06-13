%% Compute transition matrix

%  T is a 1 x A cell array. 
%  Each cell is a transisiton matrix with a numStateVector dimensional matrix.


Ttest = cell(1,sizeA);
T_gw_exp = Ttest;
for i = 1:sizeA
    T_gw_exp{i} = sparse(gw_M * exp_M, gw_M * exp_M);
    a = A(i,1);
    b = A(i,2); 
    range = gw_M*b+1:gw_M*(b+1);    % select which range is fille depending on exp state
    T_gw_exp{i}(gw_M+1:gw_M*2,gw_M+1:gw_M*2) = T_gw{a+1};  % 4th quadrant always filled
    T_gw_exp{i}(1:gw_M, range) = T_gw{a+1}; % Fill relevant 2nd range
end



for i = 1:sizeA
    tempTCell = cell(pop_M);
    for j = 1:pop_M
        for k = 1:pop_M
            tempTCell{j,k} = T_gw_exp{i} * T_pop(j,k);
        end
    end
    Ttest{i} = cell2mat(tempTCell);
end

check = Ttest{4} ~= T{4};
[rowindex, colindex] = find(check);
index = [rowindex colindex];
