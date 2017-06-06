% Inverse function of vectorIndex: takes linear state vector and returns
% index of each individual state vector

function [vectorInd] = linIndex2VecIndex(index, vectors)
    [~, numVectors] = size(vectors);
    vectorLength = cellfun(@length, vectors);
    vectorInd = zeros(numVectors,1);
    A = [1 cumprod(vectorLength)];
    indexTemp = index;
    for i = numVectors:-1:1
        vectorInd(i) = ceil(indexTemp/A(i));
        indexTemp = indexTemp - ((vectorInd(i)-1) * A(i));
    end
end