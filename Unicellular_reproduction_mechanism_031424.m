% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Spring 2024
%This simulation updates our previous simulation (Simulation 2) by not only
%visualizing clustering at a certain starvation threshold, but also having
%amoebas reproduce asexually at a prescribed iteration. This is simulated
%by having cluster sizes double at that iteration (the number of amoebas in
%a given 3.6cm grid space in the environment has doubled due to
%reproduction.)

 

%sets the size of the CA environment- how many grid squares by how many
%grid squares. Each grid square in the simulation represents 36,000
%micrometers, or, 3.6 cm based on reproduction timelines and amoeba
%movements. 
rows = 30; 
columns = 30;

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 200; %number of iterations of the simulation
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood (an 8 cell neighborhood is a 
%moore neighborhood that includes all adjacent neighbors)

reproductionTime=2; % every 10 frames amoebas will reproduce
reproductionRate = 5; %how fast amoebas reproduce(1.2 = 20% growth, 1.0 =
                        % 0% growth).


numAmoebas = 50; %number of total amoebas
amoebasPerCell = 1; %number of amoebes per grid visualization
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
 
environment =  zeros(rows, columns); % natural environment
extEnvironment = zeros(rows+2, columns+2); % environment with bounds



%create lists to update information through out simulation

environmentList = zeros(rows, columns, numIterations); %contains each cell 
        % of the environment
        
extEnvironmentList = zeros(rows+2, columns+2, numIterations);%contains the
        %environment and its boundaries
        
clusterPosList = zeros(numClusters,2,numIterations); %keeps track of the 
        %number of clusters and their positions for visualizations
        
        
clusterSizeList= zeros(numClusters,numIterations);%keeps track of all of 
        %the cluster sizes throughout the simulation (corresponds to
        %clusterPosList)

food = 2000; % Starting amount of food in the environment
starvationThreshold = 150; % Food level at which clusters start to clump
foodDecayRate = 1;  %how quickly food is being consumed on every iteration
foodList = 1:numIterations;%keeps track of how much food is available at 
        %each step in the simulation
foodList(1) = food; %sets the first value in the foodList equal to the 
                    %starting amount of food

aliveClusters = numAmoebas;%sets a value for the number of amoebas who are 
        %alive (all amoebas are alive at the beginning of the simulation)


% initialize position and size of amoeba clusters randomly
for i = 1:numClusters
    clusterPos = [randi([1 rows]) randi([1 columns])];
    environment(clusterPos(1), clusterPos(2)) = amoebasPerCell;
    clusterPosList(i,:,1) = clusterPos;
end

extEnvironment(2:rows+1, 2:columns+1) = environment;%initialize the 
        %extended environment to contain the updated environment with
        %initial clusters
        
 
% COULD BE USED TO MOVE AMOEBAS TOWARDS OTHERS
% moveTowards= @(targetPos, currentPos)(currentPos + sign(targetPos - ...
%     currentPos));


%set first index of each list correctly
environmentList(:,:,1) = environment;
extEnvironmentList(:,:,1) = extEnvironment;


%this function adds up the value of all of the neighbors in a moore
%neighborhood based on the position of a cell (x,y)
sumNeighbors = @(x, y, extEnvironment) (extEnvironment(x+1-1, y+1-1) + ...
    extEnvironment(x+1-1, y+1) + extEnvironment(x+1-1, y+1+1) + ...
    extEnvironment(x+1, y+1-1) + extEnvironment(x+1, y+1+1) + ...
    extEnvironment(x+1+1, y+1-1) + extEnvironment(x+1+1, y+1) + ...
    extEnvironment(x+1+1, y+1+1));

%Indexs for moving Amoeba clusters- keeps track of movements to access each
%neighbor in a moore neighborhood i.e. to get to the first neighbor a
%cluster would move -1 in the x direction and -1 in the y direction
indexMapping ={[-1 -1], [0 -1], [1 -1], [-1 0], [1 0], [-1 1],[0 1],[1 1]};

 
% simulation loop
for i = 2:numIterations
    %set enviorments properly for current iteration
    environment = environmentList(:,:,i-1);
    extEnvironment = extEnvironmentList(:,:,i-1);
    
    
    %cycle through each cluster and move them accordingly accordingly
    for j = 1:numClusters
        %set position and cluster size using clusterPosList which keeps
        %track of each cluster and it's current position
        clusterPos = clusterPosList(j,:,i-1);

        if(clusterPos == [0,0])%this causes the simulation to ignore 
            %clusters that have been removed from the environment after
            %combining with others
             continue;
        end
        clusterSize = environment(clusterPos(1), clusterPos(2));
         %stores the size of the cluster based on environment information
        
         
        clusterSizeList(j,i)=clusterSize;
          %adds the size of the cluster to the clusterSizeList (helps us keep
           %track of cluster sizes at each iteration of the simulation)
        
        
        environment(clusterPos(1), clusterPos(2)) = 0;
        %empties the previous location of the cluster and sets its value in
        %the environment back to 0 (there is no longer a cluster there)
        
         if rem(i,reproductionTime)==0 %check if it is a reproduction iteration
             %based on remainders from the iteration number of the
             %simulation(i) and the prescribed reproduction time
             
            clusterSize = reproductionRate *clusterSize; %increase clusterSizes 
            %based on the reproductionRate prescribed at the top
           
         end
         
         
        
        %if there are neighbors that are ameobas combine with them- use
        %sumNeighbors function to determine if neighbors values are greater
        %than 0 (empty)
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0 &&...
             food < starvationThreshold
            %get neighbors and shape into 1d array
            neighbors = reshape(extEnvironment(clusterPos(1)+1-1:...
                clusterPos(1)+1+1,clusterPos(2)-1+1:clusterPos(2)+1+1)...
                ,1,[]);
            neighbors(5) = []; % remove the clusters original position
            
            [maxNeighborSize, index] = max(neighbors);% find which direction to go
            %scan through the array of neighbors and identify the one with
            %the largest cluster size
            
            
            clusterMove = indexMapping{index};
            
            %get size of the neighboring cluster and combine the two
            environment(clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)) = maxNeighborSize + clusterSize;
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            
            %move comined cluster off the screen and remove from simulation 
             clusterPosList(j,:,i)= [0,0];
             aliveClusters = aliveClusters - 1;
        else
            
            if food < starvationThreshold && aliveClusters > 1
                % Move towards the closest cluster
                clusterMove = findClosestCluster(clusterPos, clusterPosList,...
                    i-1);
            else
                %randomly choose movement for the cluster
                indsa = randi([1, 8]);
                clusterMove = indexMapping{indsa};
            end
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows),...
                minRows), max(min(clusterPos(2) + clusterMove(2), columns), ...
                minCols)];
                    
            %update relevant variables
            environment(clusterPos(1), clusterPos(2)) = clusterSize; 
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            clusterPosList(j,:,i)= clusterPos;
            
        end
        
    end
    %update current food and food variable list
    food = cast(food - (aliveClusters * foodDecayRate),"uint8");
    foodList(i) = food;
    
    %update environments lists
    environmentList(:,:,i) = environment;
    extEnvironmentList(:,:,i) = extEnvironment;

end
 
show_CA_List(environmentList, numAmoebas, clusterPosList,clusterSizeList,...
    rows,columns,1, foodList);

%function to find closet cluster for starvation
function closestMove = findClosestCluster(clusterPos, clusterPosList, ...
    currentIteration)
    %Set min distance to largest int and direction that cluster should move
    minDistance = inf;
    closestMove = [0, 0];
    
    %go through each cluster to find the short distance between itself and
    %others
    for i = 1:size(clusterPosList, 1)
        if all(clusterPosList(i, :, currentIteration) == clusterPos) ...
                || all(clusterPosList(i, :, currentIteration) == [0, 0])
            continue; % Skip itself and clusters removed from simulation
        end
        %find distance of all other clusters
        distance = sqrt(sum((clusterPos - clusterPosList(i, :, ...
            currentIteration)).^2));
        %set min distance and direction to move if shorter distance found
        if distance < minDistance
            minDistance = distance;
            closestMove = sign(clusterPosList(i, :, currentIteration) - ...
                clusterPos);
        end
    end
end

function [ ] = show_CA_List(environmentList,numAmoebas,amoebasPosList,...
    clusterSizeList,rows,columns,interval, foodList)

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
            
            if (1-clusterSizeList(m,i)*.1)>0
                greencolor=1-clusterSizeList(m,i)*.1;
      
            else greencolor=0; 
            end
                
            
            rectangle('Position', [amoebasPosList(m,2,i)-0.5,...
                    amoebasPosList(m,1,i)-0.5, 1, 1], 'FaceColor',...
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

