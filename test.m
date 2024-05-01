environment = repmat({zeros(1,2)},5,5); % this how to create environment
% environment2 = cell(7,7);
environment2 = repmat({zeros(1,2)},7,7); % this how to create extended
                                            %environment
environmentList = cell(1,2);

%new function to get neighbor sums
sumNeighbors = @(x, y, extEnvironment) (extEnvironment{x+1-1, y+1-1}(1) + ...
    extEnvironment{x+1-1, y+1}(1) + extEnvironment{x+1-1, y+1+1}(1) + ...
    extEnvironment{x+1, y+1-1}(1) + extEnvironment{x+1, y+1+1}(1) + ...
    extEnvironment{x+1+1, y+1-1}(1) + extEnvironment{x+1+1, y+1}(1) + ...
    extEnvironment{x+1+1, y+1+1}(1));

%how to initilize environments 
for row = 1:5
   for col = 1:5
      environment{row,col} = [10 0];  
   end
end


%testing set methods
environment = setClusterSize(2,1,50, environment);
environment = setInfected(2,1,10, environment);

%setting extEnviorment with regular environment inside of it
environment2(2:6,2:6) = environment;

%testing getter methods
neighbors = getNeighborSizes(1,1,environment2);
disp(getSizeEnvironment(environment))
disp(getInfectedEnvironment(environment))

%how to set environmentlists
environmentList{1} = environment;
environmentList{2} = environment2;



function clusterSize = getClusterSize(row,col, environment)
    clusterSize = environment{row,col}(1);
end

function neighbors = getNeighborSizes(row,col,extEnviorment)
     matrixEnvironment = cell2mat(extEnviorment);
     disp(matrixEnvironment)
     neighbors = reshape(matrixEnvironment(row-1+1:...
                row+1+1,2*(col+1)-1-2:2:2*(col+1)-1+2),1,[]);
     neighbors(5) = [];
end


function clusterInfected = getInfected(row,col,environment)
    clusterInfected = environment{row,col}(2);
end

function updatedCellArray = setClusterSize(row,col, size, environment)
    environment{row,col}(1) = size;
    updatedCellArray = environment;
end

function updatedCellArray = setInfected(row,col,infected, environment)
    environment{row,col}(2) = infected;
    updatedCellArray = environment;
end

function infectedEnvironment = getInfectedEnvironment(environment)
    matrixEnvironment = cell2mat(environment);
    [numRows,numCols] = size(matrixEnvironment);
    infectedEnvironment = matrixEnvironment(:,2:2:numCols);
end

function infectedEnvironment = getSizeEnvironment(environment)
    matrixEnvironment = cell2mat(environment);
    [numRows,numCols] = size(matrixEnvironment);
    infectedEnvironment = matrixEnvironment(:,1:2:numCols-1);
end

