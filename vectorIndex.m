%% Function to get linear state index from state vector indicies

% Input
% vectors is a 1 x n cell array, where n is the number of state or action
% vectors. Each cell contains one column vectors. 
% vectorInd is a 1 x n array containing the state index for each individual
% state vector.

% Output
% index is the single linear index across the state or action space.
% S is the magnitude of the state space or action space.


function [index] = vectorIndex(vectorInd, vectors)
    [~, numVectors] = size(vectors);
    [vectorLength,~] = cellfun(@size, vectors);
    A = [1 cumprod(vectorLength)];
    index = 1:A(end);
    for i = numVectors:-1:1
        groupSize = A(i);
        groupNum = vectorInd(i);
        range = groupSize * (groupNum - 1) + 1: groupSize * groupNum;
        index = index(range);
    end

end
