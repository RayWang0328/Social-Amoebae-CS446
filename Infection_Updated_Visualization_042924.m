% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Spring 2024
%

 
rows = 30; 
columns = 30;

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 30;
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood

reproductionTime=10; % every 10 frames amoebas will reproduce
reproductionRate = 2; %how fast amoebas reproduce(1.2 = 20% growth, 1.0 =
                        % 0% growth).

   
%MAY NEED SOME OF THESE DIFFUSION CONSTANTS FOR CHEMICAL SIGNALING
% r = 0.05; % Diffusion Constant
% coolingRate = 0.2; % Cooling Constant
numAmoebas = 50; %number of total amoebas
amoebasPerCell = 5; %number of amoebes per section
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
clusterPositions = zeros(numClusters,2);

percentInfected = 0.2; %sets what proportion of amoebas are infected 0.0-1
numInfectedClusters = round(percentInfected*numClusters); %calculate number 
                                                      %of infected clusters
 
environment =  repmat({zeros(1,2)},rows,columns); % natural environment
extEnvironment = repmat({zeros(1,2)},rows+2,columns+2); % environment with bounds
%create lists to update information through out simulation
environmentList = cell(1,numIterations);
extEnvironmentList = cell(1, numIterations);
clusterPosList= cell(1,numIterations);


food = 0; % Starting amount of food in the environment
starvationThreshold = 150; % Food level at which clusters start to clump
foodDecayRate = 1;
foodList = 1:numIterations;
foodList(1) = food;

aliveClusters = numClusters;

% initialize position and size of amoeba clusters
for i = 1:numClusters
    clusterPos = [randi([1 rows]) randi([1 columns])];
    
    clusterPositions(i,:) = clusterPos;
    
    %set number of infected clusters, all amoebas in Cell will be infected
    if i<=numInfectedClusters
        environment{clusterPos(1), clusterPos(2)} = [amoebasPerCell,amoebasPerCell];
    else
        environment{clusterPos(1), clusterPos(2)} = [amoebasPerCell,0];
    end
end

%update extended enviorment
extEnvironment(2:rows+1, 2:columns+1) = environment;


%set first index of each list correctly
environmentList{1} = environment;
extEnvironmentList{1} = extEnvironment;
clusterPosList{1} = clusterPositions; 

%find the sum of the neighbors
sumNeighbors = @(x, y, extEnvironment) (extEnvironment{x+1-1, y+1-1}(1) + ...
    extEnvironment{x+1-1, y+1}(1) + extEnvironment{x+1-1, y+1+1}(1) + ...
    extEnvironment{x+1, y+1-1}(1) + extEnvironment{x+1, y+1+1}(1) + ...
    extEnvironment{x+1+1, y+1-1}(1) + extEnvironment{x+1+1, y+1}(1) + ...
    extEnvironment{x+1+1, y+1+1}(1));

%Indexs for moving Amoeba clusters
indexMapping ={[-1 -1], [0 -1], [1 -1], [-1 0], [1 0], [-1 1],[0 1],[1 1]};

 
% simulation loop
for i = 2:numIterations
    %set enviorments properly for current iteration
    environment = environmentList{i-1};
    extEnvironment = extEnvironmentList{i-1};
    oldClusterPos = clusterPosList{i-1};
    newClusterPos = [];
    
    %display for validation
%     fprintf("Iteration %d:\n",i)
%     disp(clusterCharacteristics(:,:,i-1))
    %cycle through each cluster and move them accordingly accordingly
    for j = 1:aliveClusters
        %set position and cluster size
        clusterPos = oldClusterPos(j,:);
 
        %get current characteristics of cluster
        clusterSize = getClusterSize(clusterPos(1),clusterPos(2),...
            environment);
        clusterInfected = getInfected(clusterPos(1),clusterPos(2),...
            environment);
        
        %set grid element in enviorment to be zero to show amoeabas have
        %moved
        environment{clusterPos(1),clusterPos(2)}= [0,0];
        extEnvironment(2:rows+1, 2:columns+1) = environment;
 
        
%          if rem(i,reproductionTime)==0 %check if it is a reproduction iteration
%             clusterSize = round(reproductionRate *clusterSize);
%             clusterInfected = round(reproductionRate *clusterInfected);
%          end
        
        %if there are neighbors that are ameobas, and below the starvation
        %threshold then combine with neighbors
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0 &&...
             food < starvationThreshold
         
            %get list of neighbors that is shaped into 1d array
            neighbors = getNeighborSizes(clusterPos(1),clusterPos(2),...
                extEnvironment);

            % find which direction to go and size of largest neighbor
            [maxNeighborSize, index] = max(neighbors);
            clusterMove = indexMapping{index};
            
            %Display for validation
%             fprintf("Original Cluster Size: %d, Original Cluster Inf: %d\n",...
%                 clusterSize,clusterInfected)
            

            %get size of the neighboring cluster
            neighborPos= [clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)];
         
            
            %update the clusterSize and number of infected once combined
            clusterSize = maxNeighborSize + clusterSize;
            clusterInfected =clusterInfected+getInfected(neighborPos(1),...
                neighborPos(2),environment);
            
            %print out for validation
%             fprintf("Cluster Pos: %d, %d\n",clusterPos(1),clusterPos(2))
%             fprintf("Neighbor Pos: %d, %d\n",neighborPos(1),neighborPos(2))
%             fprintf("Row Index: %d\n\n", rowIndex(1))
            

            
            %update the neighbor to show the combined size and infection
             environment = setClusterSize(neighborPos(1),...
                neighborPos(2),clusterSize,environment);
            environment = setInfected(neighborPos(1),...
                neighborPos(2),clusterInfected,environment);
%             fprintf("Updated Cluster Size: %d, UpdatedCluster Inf: %d\n\n",...
%                 clusterCharacteristics(rowIndex(1),3,i-1),...
%                 clusterCharacteristics(rowIndex(1),4,i-1))

            aliveClusters = aliveClusters - 1;
        else
            
            if food < starvationThreshold && aliveClusters > 1
                % Move towards the closest cluster
                clusterMove = findClosestCluster(clusterPos, oldClusterPos);
            else
                %randomly choose movement for the cluster
                indsa = randi([1, 8]);
                clusterMove = indexMapping{indsa};
            end
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows),...
                minRows), max(min(clusterPos(2) + clusterMove(2), columns), ...
                minCols)];          
                
            %update new position in position list
            newClusterPos= [newClusterPos; clusterPos];
            %update infection and size
            environment = setClusterSize(clusterPos(1),clusterPos(2),...
                clusterSize,environment);
            environment = setClusterSize(clusterPos(1),clusterPos(2),...
                clusterInfected,environment);
            
        end
        
    end
    %update current food and food variable list
    food = cast(food - (aliveClusters * foodDecayRate),"uint8");
    foodList(i) = food;
    
    extEnvironment(2:rows+1, 2:columns+1) = environment;
    
    %update environments lists
    environmentList{i} = environment;
    extEnvironmentList{i} = extEnvironment;
    clusterPosList{i} = newClusterPos;

end
 
show_CA_List(environmentList, numClusters,...
    rows,columns,1, foodList);

%getter and setter methods for cluster size and infection
function clusterSize = getClusterSize(row,col, environment)
    clusterSize = environment{row,col}(1);
end

function neighbors = getNeighborSizes(row,col,extEnvironment)
     matrixEnvironment = cell2mat(extEnvironment);
     neighbors = reshape(matrixEnvironment(row-1+1:...
                row+1+1,2*(col+1)-1-2:2:2*(col+1)-1+2),1,[]);
     neighbors(5) = [];
end

%function infectedEnvironment = getInfectedEnvironment(environment)
 %   matrixEnvironment = cell2mat(environment);
%    infectedEnvironment = matrixEnvironment(:,2:2:2*(columns));
%end

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

function sizeEnvironment = getSizeEnvironment(environment)
    matrixEnvironment = cell2mat(environment);
    [numRows,numCols] = size(matrixEnvironment);
    sizeEnvironment = matrixEnvironment(:,1:2:numCols-1);
end

%function to find closet cluster for starvation
function closestMove = findClosestCluster(clusterPos, clusterPosList)
    %Set min distance to largest int and direction that cluster should move
    minDistance = inf;
    closestMove = [0, 0];
    
    %go through each cluster to find the short distance between itself and
    %others
    for i = 1:size(clusterPosList, 1)
        if all(clusterPosList(i, :) == clusterPos) ...
                || all(clusterPosList(i, :) == [0, 0])
            continue; % Skip itself and clusters removed from simulation
        end
        %find distance of all other clusters
        distance = sqrt(sum((clusterPos - clusterPosList(i, :)).^2));
        %set min distance and direction to move if shorter distance found
        if distance < minDistance
            minDistance = distance;
            closestMove = sign(clusterPosList(i, :) - ...
                clusterPos);
        end
    end
end

%don't think you need cluster characteristics should just need the environment
% take enviroment list and cycle through it using 
function [ ] = show_CA_List(environmentList,numAmoebas,...
    rows,columns,interval, foodList)

    for i=1:interval:length(environmentList)
        environment = environmentList{i};
        infectedEnvironment=getInfectedEnvironment(environment)
        sizeEnvironment=getSizeEnvironment(environment)
        

       %map to set colors for legend
        redMap=[ 1 1 1
             .1,0,0
             .15,0,0
             .2,0,0
             .25,0,0
             .3,0,0
             .35,0,0
             .4,0,0
             .45,0,0
             .5,0,0
             .55,0,0
             .6,0,0
             .65,0,0
             .7,0,0
             .75,0,0
             .8,0,0
             .85,0,0
             .9,0,0];
         
         
         greenMap=[ 1 1 1
              0,.1,0
              0,.15,0
              0,.2,0
              0,.25,0
              0,.3,0
              0,.35,0
              0,.4,0
              0,.45,0
              0,.5,0
              0,.55,0
              0,.5,0
              0,.65,0
              0,.7,0
              0,.75,0
              0,.8,0
              0,.85,0
              0,.9,0];
         
 
        
        subplot(1,2,1); 
        colormap(redMap); 
        imagesc(infectedEnvironment);
        caxis([0,18]);
        colorbar;
        bar1=colorbar;
        bar1.Ticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];   
        bar1.TickLabels={'empty','1 infected amoeba','2 infected amoeba ',...
            '3 infected amoeba','4 infected amoeba','5 infected amoeba',...
            '6 infected amoeba','7 infected amoeba','8 infected amoeba',...
            '9 infected amoeba','10 infected amoeba','11 infected amoeba',...
            '12 infected amoeba','13 infected amoeba','14 infected amoeba',...
            '15 infected amoeba','16 infected amoeba','17 infected amoeba',...
            '18+ infected amoeba'};
        
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Infected Environment - Frame: %d, Food Available: %d', i, foodList(i)));
        set(gca,'YDir','reverse');

        subplot(1,2,2); 
        colormap(greenMap); 
        imagesc(sizeEnvironment);
        caxis([0,12]);
        colorbar;
        bar2=colorbar;
        bar2.Ticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];   
        bar2.TickLabels={'empty','1 amoeba','2 amoeba ',...
            '3 amoeba','4 amoeba','5 amoeba',...
            '6 amoeba','7 amoeba','8 amoeba',...
            '9 amoeba','10 amoeba','11 amoeba',...
            '12 amoeba','13 amoeba','14 amoeba',...
            '15 amoeba','16 amoeba','17 amoeba',...
            '18+ amoeba'};
        
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Size Environment - Frame: %d, Food Available: %d', i, foodList(i)));
        set(gca,'YDir','reverse');

        w = waitforbuttonpress;
    end
end