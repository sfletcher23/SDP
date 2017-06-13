function [T_gw_new] = gw_transmat_multiwell(T_gw, gw_M, gw_M_well)

T_gw_new = cell(1,1);

% A = [0 0]
T_gw_new{1} = spdiags(ones(gw_M,1),0,gw_M,gw_M);

% A = [1 0]
T_gw_new{2} = sparse(gw_M,gw_M);
for i = 1:gw_M_well
    range = (i-1)*gw_M_well+1: i*gw_M_well;
    T_gw_new{2}(range,range) = T_gw{2};
end

% A = [0 1]
T_gw_new{3} = sparse(gw_M,gw_M);
for i = 1:gw_M_well
    range1 = i:gw_M_well:gw_M_well * (gw_M_well-1)+i;
    range2 = gw_M_well+i:gw_M_well: gw_M;
    T_gw_new{3}(range1,range1) = T_gw{2};
end

% A = [1 1]
T_gw_new{4} = sparse(gw_M,gw_M);
temp1 = sparse(gw_M,gw_M);
temp2 = sparse(gw_M,gw_M);
for i = 1:gw_M_well
    range1 = (i-1)*gw_M_well+1: i*gw_M_well;
    if i<gw_M_well
        range2 = range1 + gw_M_well;
    else
        range2 = range1;
    end
   temp1(range1,range1) = T_gw{2} ./2;
   temp2(range1,range2) = T_gw{2} ./2;
   T_gw_new{4} = temp1 + temp2;
end

% Check validity of transition matrices
for i = 1:length(T_gw_new)
    rowsum = sum(T_gw_new{i},2);
    err = abs(rowsum - 1);
    if err > 0
        error(['Invalid transition matrix action ', num2str(i)])
    end
end