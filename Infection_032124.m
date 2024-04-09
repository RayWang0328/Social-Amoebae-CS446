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
numIterations = 15;
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

percentInfected = 0.2; %sets what proportion of amoebas are infected 0-1
%calculate number of infected clusters.
numInfectedClusters = round(percentInfected*numClusters);
 
environment =  zeros(rows, columns); % natural environment
extEnvironment = zeros(rows+2, columns+2); % environment with bounds

%create lists to update information through out simulation
environmentList = zeros(rows, columns, numIterations);
extEnvironmentList = zeros(rows+2, columns+2, numIterations);


%first two numbers are x,y coordinates 
%third number represents total cluster size and fourth number is the amount
%of infected amoeba in the cluster
clusterCharacteristics= zeros(numClusters,4, numIterations);

food = 0; % Starting amount of food in the environment
starvationThreshold = 150; % Food level at which clusters start to clump
foodDecayRate = 1;
foodList = 1:numIterations;
foodList(1) = food;

aliveClusters = numClusters;

% initialize position and size of amoeba clusters
for i = 1:numClusters
    clusterPos = [randi([1 rows]) randi([1 columns])];
    environment(clusterPos(1), clusterPos(2)) = amoebasPerCell;
    clusterCharacteristics(i,1:2,1) = clusterPos;
    
    %set number of infected clusters
    if i<=numInfectedClusters
        clusterCharacteristics(i,3:4,1) = [amoebasPerCell,amoebasPerCell];
    else
        clusterCharacteristics(i,3:4,1) = [amoebasPerCell,0];
    end
end

extEnvironment(2:rows+1, 2:columns+1) = environment;


%set first index of each list correctly
environmentList(:,:,1) = environment;
extEnvironmentList(:,:,1) = extEnvironment;

%find the sum of the neighbors
sumNeighbors = @(x, y, extEnvironment) (extEnvironment(x+1-1, y+1-1) + ...
    extEnvironment(x+1-1, y+1) + extEnvironment(x+1-1, y+1+1) + ...
    extEnvironment(x+1, y+1-1) + extEnvironment(x+1, y+1+1) + ...
    extEnvironment(x+1+1, y+1-1) + extEnvironment(x+1+1, y+1) + ...
    extEnvironment(x+1+1, y+1+1));

%Indexs for moving Amoeba clusters
indexMapping ={[-1 -1], [0 -1], [1 -1], [-1 0], [1 0], [-1 1],[0 1],[1 1]};

 
% simulation loop
for i = 2:numIterations
    %set enviorments properly for current iteration
    environment = environmentList(:,:,i-1);
    extEnvironment = extEnvironmentList(:,:,i-1);
    clusterPositions = clusterCharacteristics(:,1:2,i-1);
    
    fprintf("Iteration %d:\n",i)
    disp(clusterCharacteristics(:,:,i-1))
    %cycle through each cluster and move them accordingly accordingly
    for j = 1:numClusters
        %set position and cluster size
        clusterPos = clusterPositions(j,:);

        if(clusterPos == [0,0])
             continue;
        end
        clusterSize = clusterCharacteristics(j,3,i-1);
        clusterInfected = clusterCharacteristics(j,4,i-1);

        environment(clusterPos(1), clusterPos(2)) = 0;
        extEnvironment(clusterPos(1), clusterPos(2)) = 0;
        neighbors = reshape(extEnvironment(clusterPos(1)+1-1:...
                clusterPos(1)+1+1,clusterPos(2)-1+1:clusterPos(2)+1+1)...
                ,1,[]);
        
%          if rem(i,reproductionTime)==0 %check if it is a reproduction iteration
%             clusterSize = round(reproductionRate *clusterSize);
%             clusterInfected = round(reproductionRate *clusterInfected);
%          end
        
        %if there are neighbors that are ameobas combine with them
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0 &&...
             food < starvationThreshold
            %get neighbors and shape into 1d array
            neighbors = reshape(extEnvironment(clusterPos(1)+1-1:...
                clusterPos(1)+1+1,clusterPos(2)-1+1:clusterPos(2)+1+1)...
                ,1,[]);
            neighbors(5) = []; % remove the clusters original position
            [maxNeighborSize, index] = max(neighbors);% find which direction to go
            clusterMove = indexMapping{index};
            
            %update cluster size from max neighbor
            clusterSize = maxNeighborSize + clusterSize;
            
            
        
            %get size of the neighboring cluster and combine the two
            neighborPos= [clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)];
            

            %find number of infected
            rowIndex = find(ismember(clusterPositions,neighborPos, 'rows'));   
            clusterInfected = clusterInfected + clusterCharacteristics...
                (rowIndex,4,i-1);

            fprintf("Cluster Pos: %d, %d\n",clusterPos(1),clusterPos(2))
            fprintf("Neighbor Pos: %d, %d\n",neighborPos(1),neighborPos(2))
            fprintf("Row Index: %d\n\n", rowIndex)
            
            
            %move comined cluster off the screen and remove from simulation 
            clusterCharacteristics(j,1:2,i)= [0,0];
            clusterPositions(j,:) = [0,0];
            
            fprintf("Cluster Size: %d, Cluster Inf: %d\n",clusterSize,clusterInfected)
            clusterCharacteristics(rowIndex,3:4,i)=[clusterSize,clusterInfected];


            aliveClusters = aliveClusters - 1;
        else
            
            if food < starvationThreshold && aliveClusters > 1
                % Move towards the closest cluster
                clusterMove = findClosestCluster(clusterPos, clusterPositions);
            else
                %randomly choose movement for the cluster
                indsa = randi([1, 8]);
                clusterMove = indexMapping{indsa};
            end
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows),...
                minRows), max(min(clusterPos(2) + clusterMove(2), columns), ...
                minCols)];          
                
            
            clusterCharacteristics(j,1:2,i)= clusterPos;
            %update infection and size
            clusterCharacteristics(j,3,i) = clusterSize;
            clusterCharacteristics(j,4,i) = clusterInfected;
            
        end
        
    end
    %update current food and food variable list
    food = cast(food - (aliveClusters * foodDecayRate),"uint8");
    foodList(i) = food;
    
    %update enviorment with new cluster sizes
    for j = 1:numClusters
        if(clusterCharacteristics(j,1,i) ~= 0)
        environment(clusterCharacteristics(j,1,i),...
            clusterCharacteristics(j,2,i)) = clusterCharacteristics(j,3,i);
        end
    end
    extEnvironment(2:rows+1, 2:columns+1) = environment;
    
    %update environments lists
    environmentList(:,:,i) = environment;
    extEnvironmentList(:,:,i) = extEnvironment;

end
 
show_CA_List(environmentList, numClusters, clusterCharacteristics,...
    rows,columns,1, foodList);

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

function [ ] = show_CA_List(environmentList,numAmoebas,clusterCharacteristics,...
    rows,columns,interval, foodList)

    for i=1:interval:length(environmentList)
        environment = environmentList(:,:,i);
        hold on;

       %map to set colors for legend
        map=[ 0.478 0.318 0.102
             0.478 0.318 0.102
             0,.9,0
             0,.8,0
             0,.7,0
             0,.6,0
             0,.5,0
             0,.4,0
             0,.3,0
             0,.2,0
             0,.1,0
             0,0,0];
        imagesc(environment);
        colormap(map);
        
        %more set up for legend
        caxis([0,12]);
        lifeCycleColors=colorbar;
        lifeCycleColors.Ticks=[1,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5];   
        lifeCycleColors.TickLabels={'empty','1 amoeba','2 amoeba cluster',...
            '3 amoeba cluster','4 amoeba cluster','5 amoeba cluster',...
            '6 amoeba cluster','7 amoeba cluster','8 amoeba cluster',...
            '9 amoeba cluster','10+ amoeba cluster'};
         hold;
        
        
       
        % Plot mammals positions on top of the heat map
        for m = 1:numAmoebas
            
            if (1-clusterCharacteristics(m,3,i)*.1)>0
                greencolor=1-clusterCharacteristics(m,3,i)*.1;
      
            else greencolor=0; 
            end
                
            
            rectangle('Position', [clusterCharacteristics(m,2,i)-0.5,...
                    clusterCharacteristics(m,1,i)-0.5, 1, 1], 'FaceColor',...
                    [0,greencolor,0],'EdgeColor', 'k');
                
        end    
           
             


        hold off;
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Simulation Frame: %d, Food Available: %d', i, ...
            foodList(i)));
        set(gca,'YDir','reverse');

        w=waitforbuttonpress;
        
    end
end