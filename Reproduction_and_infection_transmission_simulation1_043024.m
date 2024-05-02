% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Spring 2024


% This is an implementation of the model of Dictyostelium discoideum amoeba social
% behavior. Social amoeba behavior depends on the condition of starvation
% in the simulation. When a starvation threshold is reached, clustering
% behavior occurs which causes amoeba to clump together into a slug. 
% Reproduction occurs for unicellular amoeba before completely aggregated
% into a multicellular slug. Reproduction rates can be set to different
% values for infected and uninfected amoebas as per *biological evidence*
% that infection may impact reproduction rates. Amoeba
% can be infected by the bacteria Burkholderia at prescribed rates.
% Noh et al's "Facultative symbiont virulence determines horizontal ...
% transmission rate without host specificity in Dictyostelium discoideum ...
% social amoebas"  indicates that infection spread has a linear
% relationship with the amount of infected amoeba in the local environment.
% As such, increase in infected amoeba is modeled as proportional to the
% density of infection in the surrounding "cluster." As per *biological
% evidence* around 25-35% of infected amoeba persist in the population.
% Given this we can assume that around 65% of amoeba eventually die off
% from their infection. This is modeled with an infectionDeathRate that
% occurs after a certain number of hours of being infected. Further social 
% amoeba behavior involves a spore and fruiting body dispersal mechanism 
% following slug formation. This is implemented as well. Random seeds are 
% used to store data for analysis. Potential
% future uses for this implementation include running simulations with
% varying reproduction and infection rates in order to understand what
% causes an equilibrium with persistent infection, as well as explore the
% effect of reproductive rates of both the infected and uninfected
% populations on the spread of infection during social clustering behavior.


%Current biological understanding has amoeba reproducing every 4 hours as 
%per Professor Suegene Noh. Amoeba can move around 10 micrometers per
%second. If we use each frame of a simulation to represent 1 hour, and
%amoeba move one grid square on each frame, then each grid square
%represents 36,000 micrometers or 3.6 centimeters. 

%these values represent the dimensions of the 2D environment we are
%simulating. 
rows = 30; %the width of the environment- each row increases width by 3.6 ...
                                                                        %cm
columns = 30; % the height of the environment- each row increases height...
                                                                %by 3.6 cm

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables:

%This is the length that the simulation will run for. Each iteration...
%represents an hour of time passing in the simulation
numIterations = 100;


x = 1:rows;
y = 1:columns;

%This is a parameter that can be used to turn integer visualization on and
%off. A value of 0 will visualize environments with only colored squares. A
%value of 1 will visualize environments with overlayed integer values for
%infected number and cluster size- for better readability. 
integerVisualization=1;


neighborhoodSize = 8; % size of neighborhood. In the cellular automata...
%model, moore neighborhoods mean that the 8 adjacent cells will be
%registered when updating the current cell. 

%this sets the interval at which amoebas will reproduce. Standard
%biological value would indicate that amoebas reproduce every 4 hours, so
%that value is set here. This can be changed depending on environmental
%conditions. 
reproductionTime=4; 

%This value determines the interval at which infection will spread in the
%population. Given that Burkholderia can reproduce around every 20 minutes,
%we anticipate that infection would spread horizontally  in adjacent 
%(clustered) amoebas on every iteration of the simulation. (Every 1 hour)
infectionTime=1;%every 2 hours infection will spread in the clusters

%At a certain point, infected amoeba will succumb to their infection and
%die. This value indicates how frequently in the simulation that happens.
%Can be changed by domain experts to increase biological accuracy. 
infectionLimit=12; %after 12 hours of persistent infection, amoeba will die

%There is biological evidence per Professor Suegene Noh that different
%reproduction rates may exist for infected and uninfected amoeba. This is
%implemented here. These values can be changed by domain experts to
%increase biological accuracy. 
infectedReproductionRate = 2; %how fast amoebas reproduce(1.2 = 20% growth, 1.0 =
uninfectedReproductionRate=1.2;   % 0% growth).

%This indicates the percentage of amoebas that will die off when we reach
%the infection limit. 25-35% of amoebas are believed to persist with
%infection so we assume that around 65% of them will succumb to their
%infection and die. This value can be changed by the domain expert to
%increase biological accuracy. 
infectionDeathRate=.65;

%this is the starting number of amoebas in the simulation.  
numAmoebas = 1000; %number of initialtotal amoebas

%Given the magnitude of amoebas in a simulation, agent based modeling for
%each individual amoeba would be computationally intensive. As such,
%aggregate modeling was chosen in order to better visualize large numbers
%of these unicellular organism. 
amoebasPerCell = 10; %number of amoebes per grid square (how we want to 
%initially divide up our amoeba population into aggregate "clusters"

numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
clusterPositions = zeros(numClusters,2); %sets a storage location for 
                                            %cluster positions

percentInfected = 0.2; %sets what proportion of amoebas are infected 0.0-1
numInfectedClusters = round(percentInfected*numClusters); %calculate number 
                                                      %of infected clusters
                                                      
                                                      
numInfected=percentInfected*numAmoebas; %use the initial percent infected 
%value to calculated the number of infected amoeba at the beginning of the 
%simulation. 


 
environment =  repmat({zeros(1,2)},rows,columns); % natural environment
%this consists of a cell array containing cells with 2 values in them.
%These will later be updated to store infected amoebas and cluster sizes


extEnvironment = repmat({zeros(1,2)},rows+2,columns+2); % environment with bounds
%part of the cellular automata process to ensure edge cases can read in
%neighbors and update. 

%create lists to update information through out simulation. this is used
%for visualization of each iteration of the simulation. 
environmentList = cell(1,numIterations); %stores each environment
extEnvironmentList = cell(1, numIterations); %stores each extended 
%environment
clusterPosList= cell(1,numIterations); %stores all cluster position lists


%Amoeba social behavior is dependent on starvation signalling which
%triggers amoebas to cluster together. As such, there is  prescribed amount
%of food in the environment at the beginning of the simulation which is
%steadily consumed by the amoeba in the environment. When this food
%decreases to a certain value or below a "starvation threshold" and can no
%longer support the population of amoebas, clustering behavior is
%triggered.
food = 100; % Starting amount of food in the environment
starvationThreshold = 150; % Food level at which clusters start to clump
foodDecayRate = 1;%rate at which amoebas consume available food
foodList = 1:numIterations; %stores food availability at each point in 
                                                        %the simulation
foodList(1) = food; %sets the first stored food value to the initial allotted
%food for the environment

aliveClusters = numClusters; %keeps track of the number of clusters that are 
%active in the environment. All clusters are initially active. 

% initialize position and size of amoeba clusters
for i = 1:numClusters
    
    %sets the clusters position randomly in the envrionment
    clusterPos = [randi([1 rows]) randi([1 columns])];
    
    %stores randomly generated cluster positions
    clusterPositions(i,:) = clusterPos;
    
    %set number of infected clusters, all amoebas in cluster will be infected
    if i<=numInfectedClusters
        
        %sets clusterInfected to all amoebas in the cluster 
        environment{clusterPos(1), clusterPos(2)} = [amoebasPerCell,amoebasPerCell];
    else
        %initializes the cluster with no infected amoeba
        environment{clusterPos(1), clusterPos(2)} = [amoebasPerCell,0];
    end
end

%update extended environment to contain the newly initiated environment
extEnvironment(2:rows+1, 2:columns+1) = environment;


%set first index of each storage list correctly
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
    
   
    %cycle through each cluster and move them accordingly accordingly
    for j = 1:aliveClusters
        %set position and cluster size
        clusterPos = oldClusterPos(j,:);
 
        %get current characteristics of cluster
        clusterSize = getClusterSize(clusterPos(1),clusterPos(2),...
            environment);
        clusterInfected = getInfected(clusterPos(1),clusterPos(2),...
            environment);
        clusterUninfected=clusterSize-clusterInfected;
        
        %set grid element in environment to be zero to show amoeabas have
        %moved
        environment{clusterPos(1),clusterPos(2)}= [0,0];
        extEnvironment(2:rows+1, 2:columns+1) = environment;
 
        
        %reproduce on reproductionTime
        if rem(i,reproductionTime)==0 && aliveClusters>1 %check if it is a reproduction iteration
            
            
            clusterInfected = round(infectedReproductionRate *clusterInfected);
             
            clusterUninfected=round(uninfectedReproductionRate*clusterUninfected);
            
            clusterSize = clusterInfected+clusterUninfected;
            
            
           
             
        end
        
        %transmit infection on infectiontime
        if rem(i,infectionTime)==0 && clusterSize~=0
           infectedPercentage=(clusterInfected/clusterSize); 
           clusterInfected=clusterInfected + round(clusterUninfected*infectedPercentage);
            
           clusterInfected=round(clusterInfected);
           
           
           
           
        end
        
        
        %25-35% of infected amoebas persist with infection so around .65%
        %will die when they reach the infection limit iteration
        
        if rem(i,infectionLimit)==0
          infectionKilled=infectionDeathRate*clusterInfected;
          clusterInfected=clusterInfected-infectionKilled;
          clusterSize=clusterSize-infectionKilled;
          
          clusterInfected=round(clusterInfected);
          clusterSize=round(clusterSize);
        
        end
        
        
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
            
            

            %get size of the neighboring cluster
            neighborPos= [clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)];
         
            
            %update the clusterSize and number of infected once combined
            clusterSize = maxNeighborSize + clusterSize;
            clusterInfected =clusterInfected+getInfected(neighborPos(1),...
                neighborPos(2),environment);
            
           

            
            %update the neighbor to show the combined size and infection
             environment = setClusterSize(neighborPos(1),...
                neighborPos(2),clusterSize,environment);
            environment = setInfected(neighborPos(1),...
                neighborPos(2),clusterInfected,environment);
%             

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
            environment = setInfected(clusterPos(1),clusterPos(2),...
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
 
show_CA_List(environmentList, ...
    rows,columns,1, foodList,integerVisualization,numAmoebas,numInfected);

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


function [ ] = show_CA_List(environmentList,...
    rows,columns,interval, foodList,integerVisualization,numAmoebas,numInfected)
    
    infectedInterval=numInfected/20;
    sizeInterval=numAmoebas/20;
    
    for i=1:interval:length(environmentList)
        environment = environmentList{i};
        infectedEnvironment=getInfectedEnvironment(environment);
        sizeEnvironment=getSizeEnvironment(environment);
        hold on;
        
        
        subplot(1,2,1); 
        colormap(subplot(1,2,1),flipud(bone)); 
        imagesc(infectedEnvironment);
        caxis([0,numInfected]);
        colorbar;
        bar1=colorbar;
        bar1.Ticks=[infectedInterval,2*infectedInterval,3*infectedInterval...
            ,4*infectedInterval,5*infectedInterval,6*infectedInterval,...
            7*infectedInterval,8*infectedInterval,9*infectedInterval,...
            10*infectedInterval,11*infectedInterval,12*infectedInterval...
            ,13*infectedInterval,14*infectedInterval,15*infectedInterval,...
            16*infectedInterval,17*infectedInterval,18*infectedInterval,...
            19*infectedInterval,20*infectedInterval,21*infectedInterval]; 
        
        bar1.TickLabels={'empty',sprintf('%d infected amoeba',infectedInterval),...
            sprintf('%d infected amoeba',2*infectedInterval),...
            sprintf('%d infected amoeba',3*infectedInterval),...
            sprintf('%d infected amoeba',4*infectedInterval),...
            sprintf('%d infected amoeba',5*infectedInterval),...
            sprintf('%d infected amoeba',6*infectedInterval),...
            sprintf('%d infected amoeba',7*infectedInterval),...
            sprintf('%d infected amoeba',8*infectedInterval),...
            sprintf('%d infected amoeba',9*infectedInterval),...
            sprintf('%d infected amoeba',10*infectedInterval),...
            sprintf('%d infected amoeba',11*infectedInterval),...
            sprintf('%d infected amoeba',12*infectedInterval),...
            sprintf('%d infected amoeba',13*infectedInterval),...
            sprintf('%d infected amoeba',14*infectedInterval),...
            sprintf('%d infected amoeba',15*infectedInterval),...
            sprintf('%d infected amoeba',16*infectedInterval),...
            sprintf('%d infected amoeba',17*infectedInterval),...
            sprintf('%d infected amoeba',18*infectedInterval),...
            sprintf('%d infected amoeba',19*infectedInterval),...
            sprintf('%d infected amoeba',20*infectedInterval)};
        
           for m = 1:size(infectedEnvironment)[1];
                for j= 1:size(infectedEnvironment)[2];

            
                    if (infectedEnvironment(m,j)==0 && sizeEnvironment(m,j)>0)
                        subplot(1,2,1);
                        rectangle('Position', [j-0.5,m-0.5, 1, 1], 'FaceColor',...
                        [1,1,1],'EdgeColor', 'k');
                    
                        text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'red');
            
                    end
                end       
           end  
           
        for m = 1:size(infectedEnvironment, 1)
            for j = 1:size(infectedEnvironment, 2)
                    if (integerVisualization==1 && infectedEnvironment(m,j)~=0)
                    text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'red');
                    end
            end
        end
        
        hold off;
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Infected Environment - Frame: %d, Food Available: %d', i, foodList(i)));
        set(gca,'YDir','reverse');

        
        hold on;
        colormap(subplot(1,2,2),parula);
        subplot(1,2,2);  
        imagesc(sizeEnvironment);
        caxis([0,numAmoebas]);
        colorbar;
        bar2=colorbar;
        
        
        bar2.Ticks=[sizeInterval,2*sizeInterval,3*sizeInterval...
            ,4*sizeInterval,5*sizeInterval,6*sizeInterval,...
            7*sizeInterval,8*sizeInterval,9*sizeInterval,...
            10*sizeInterval,11*sizeInterval,12*sizeInterval...
            ,13*sizeInterval,14*sizeInterval,15*sizeInterval,...
            16*sizeInterval,17*sizeInterval,18*sizeInterval,...
            19*sizeInterval,20*sizeInterval,21*sizeInterval]; 
        
        
        bar2.TickLabels={'empty',sprintf('%d infected amoeba',infectedInterval),...
            sprintf('%d infected amoeba',2*sizeInterval),...
            sprintf('%d infected amoeba',3*sizeInterval),...
            sprintf('%d infected amoeba',4*sizeInterval),...
            sprintf('%d infected amoeba',5*sizeInterval),...
            sprintf('%d infected amoeba',6*sizeInterval),...
            sprintf('%d infected amoeba',7*sizeInterval),...
            sprintf('%d infected amoeba',8*sizeInterval),...
            sprintf('%d infected amoeba',9*sizeInterval),...
            sprintf('%d infected amoeba',10*sizeInterval),...
            sprintf('%d infected amoeba',11*sizeInterval),...
            sprintf('%d infected amoeba',12*sizeInterval),...
            sprintf('%d infected amoeba',13*sizeInterval),...
            sprintf('%d infected amoeba',14*sizeInterval),...
            sprintf('%d infected amoeba',15*sizeInterval),...
            sprintf('%d infected amoeba',16*sizeInterval),...
            sprintf('%d infected amoeba',17*sizeInterval),...
            sprintf('%d infected amoeba',18*sizeInterval),...
            sprintf('%d infected amoeba',19*sizeInterval),...
            sprintf('%d + infected amoeba',20*sizeInterval)};
        
       
       for m = 1:size(sizeEnvironment, 1)
            for j = 1:size(sizeEnvironment, 2)
                if (integerVisualization==1 && sizeEnvironment(m,j)~=0)
                    text(j, m, num2str(sizeEnvironment(m, j)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'k');
                end
            end
        end
        
        
        
        hold off;
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Size Environment - Frame: %d, Food Available: %d', i, foodList(i)));
        set(gca,'YDir','reverse');

        w = waitforbuttonpress;
    end
end