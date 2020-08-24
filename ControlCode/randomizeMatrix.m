% function to generate a matrix by randomizing columns in input matrix
function randMatrix = randomizeMatrix(inputMatrix)
    nRows = size(inputMatrix, 1);
    nColumns = size(inputMatrix, 2);
    randMatrix = zeros(nRows, nColumns);
    for i = 1:nRows
        random_idx = randperm(nColumns);
        randMatrix(i, :) = inputMatrix(i, random_idx);
    end
end