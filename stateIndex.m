%% Function to get linear state index from state vector indicies

% Input
% vectors is a 1 x n cell array, where n is the number of state or action
% vectors. Each cell contains one column vectors. 
% vectorInd is a 1 x n array containing the state index for each individual
% state vector.

% Output
% index is the single linear index across the state or action space.
% N is the magnitude of the state space or action space.


function [index, S] = stateIndex(vectorInd, vectors)
    [~, numVectors] = size(vectors);
    S = 1;
    for i = 1:numVectors
        S = S * length(vectors{i});
    end
    S_size = [];
    for i = 1:numVectors
        S_size(end+1) = S;
    end
    vectorIndCell = num2cell(vectorInd);
    index = sub2ind(S_size,vectorIndCell{:});

end
