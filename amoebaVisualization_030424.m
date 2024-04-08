% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Fall, 2023
% This is an update to simulation 1 where clustering occurs based on random
% amoeba movements. Clusters are visualized using increasing green color
% based on amoeba concentrations. 

 %sets the size of the CA environment- how many grid squares by how many
%grid squares. Each grid square in the simulation represents 36,000
%micrometers, or, 3.6 cm based on reproduction timelines and amoeba
%movements.  
rows = 20; 
columns = 20;

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 200;%number of iterations the simulation will go
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood (an 8 cell neighborhood is a 
%moore neighborhood that includes all adjacent neighbors)


numAmoebas = 50; %number of total amoebas in the simulation
amoebasPerCell = 1; %number of amoebes per grid visualization
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
 
environment =  zeros(rows, columns); % natural environment
extEnvironment = zeros(rows+2, columns+2); % environment with bounds



%create lists to update information through out simulation


environmentList = zeros(rows, columns, numIterations); %contains each grid 
        %of the environment

extEnvironmentList = zeros(rows+2, columns+2, numIterations);%contains the
        %environment and its boundaries

clusterPosList = zeros(numClusters,2,numIterations);  %keeps track of the 
        %number of clusters and their positions for visualization

clusterSizeList= zeros(numClusters,numIterations); %keeps track of all of 
        %the cluster sizes throughout the simulation (corresponds to
        %clusterPosList)




% initialize position and size of amoeba clusters randomly
for i = 1:numClusters
    clusterPos = [randi([1 rows]) randi([1 columns])];
    environment(clusterPos(1), clusterPos(2)) = amoebasPerCell;
    clusterPosList(i,:,1) = clusterPos;
end

extEnvironment(2:rows+1, 2:columns+1) = environment; %initialize the 
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
    
 
   
    %cycle through each cluster and move them accordingly 
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
        
        %if there are neighbors that are ameobas combine with them- use
        %sumNeighbors function to determine if neighbors values are greater
        %than 0 (empty)
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0
            
            %get neighbors and shape into 1d array
            neighbors = reshape(extEnvironment(clusterPos(1)+1-1:...
                clusterPos(1)+1+1,clusterPos(2)-1+1:clusterPos(2)+1+1)...
                ,1,[]);
            neighbors(5) = []; % remove the clusters original position
            
            [maxNeighborSize, index] = max(neighbors);% find which direction to go
            %scan through the array of neighbors and identify the one with
            %the largest cluster size
            
            clusterMove = indexMapping{index};
            %based on identifiying cluster size and neighbor position, use
            %the index mapping array to direct the cluster which direction
            %to move based on the index of the largest neighbor
            
            
            %get size of the neighboring cluster and combine the two, this
            %updates the environment to reflect the size of the amoeba
            %clusters
            environment(clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)) = maxNeighborSize + clusterSize;
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            
            %move comined cluster off the screen and remove from simulation
             %This effectively visualizes how clusters have combined into
            %one. 
             clusterPosList(j,:,i)= [0,0];
        else
            %randomly choose movement for the cluster if it has no amoeba
            %neighbors (all neighbors are empty)
            
            %determine random index for movement directions
            indsa= randi([1, 8]);
            
            %assign random movement using indexMapping array
            clusterMove = indexMapping{indsa};
            
            %update cluster position to reflect this movement
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows), minRows), ...
                        max(min(clusterPos(2) + clusterMove(2), columns), minCols)];
                    
            %update relevant variables
            
            %change environment cell in grid to reflect the size of the
            %cluster it holds
            environment(clusterPos(1), clusterPos(2)) = clusterSize; 
            
            %update the extended environment to contain the environment
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            
            
            %update the cluster position list to contain the new cluster
            %position after movement
            clusterPosList(j,:,i)= clusterPos;
            
        end
        
    end
    
    
    %add the updated environment to the environmentList
    environmentList(:,:,i) = environment;
    
    %add the updated extended environment to the extended environment list.
    extEnvironmentList(:,:,i) = extEnvironment;

end
 

%feeds the environment List, number of Amoebas, cluster position list,
%rows, columns and interval to the function show_CA_list for visualization
%of our simulation
show_CA_List(environmentList, numAmoebas, clusterPosList,clusterSizeList,rows,columns,1);



%a functon for visualizing each frame in the simulation> 
function [ ] = show_CA_List(environmentList,numAmoebas,amoebasPosList,clusterSizeList,rows,columns,interval)

    for i=1:interval:length(environmentList)
        %for each frame in the simulation, determine the environment at
        %that time step
        environment = environmentList(:,:,i);
        hold on;

       %map to set colors for legend. This visualizes increasing density of
       %amoeba clusters where 0 is an empty cell and increasing
       %clusterSizes have increasing shades of green
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
        
        %More set up for legend
        caxis([0,12]);
        lifeCycleColors=colorbar;
        
        %set ticks for better colorbar visualization
        lifeCycleColors.Ticks=[1,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5];  
        
        %label ticks to reflect the sizes of visualized clusters
        lifeCycleColors.TickLabels={'empty','1 amoeba','2 amoeba cluster','3 amoeba cluster',...
        '4 amoeba cluster','5 amoeba cluster','6 amoeba cluster','7 amoeba cluster','8 amoeba cluster','9 amoeba cluster','10+ amoeba cluster'};
         hold;
        
        
       
        % Plot clusters positions on top of the environment
        for m = 1:numAmoebas
            
            %if the cluster at this point contains amoebas (is not 0)
            if (1-clusterSizeList(m,i)*.1)>0
                greencolor=1-clusterSizeList(m,i)*.1;%determine the 
                %green concentration to be visualized 
      
            else greencolor=0; %if no amoebas are available at this cluster,
                %visualize the cluster with 0 green color
            end
                
            
            rectangle('Position', [amoebasPosList(m,2,i)-0.5,...
                    amoebasPosList(m,1,i)-0.5, 1, 1], 'FaceColor',[0,greencolor,0], ...
                    'EdgeColor', 'k');
                
        end    
           
             


        hold off;
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Simulation Frame: %d', i));
        set(gca,'YDir','reverse');

        w=waitforbuttonpress;
        
    end
end

