% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Spring 2024


% This is an implementation of the model of Dictyostelium discoideum amoeba
% social behavior.
% Social amoeba behavior depends on the condition of starvation
% in the simulation. When a starvation threshold is reached, clustering
% behavior occurs which causes amoeba to clump together into a slug. 
% Reproduction occurs for unicellular amoeba before completely aggregated
% into a multicellular slug. Reproduction rates can be set to different
% values for infected and uninfected amoebas as per *biological evidence*
% that infection may impact reproduction rates. Amoeba
% can be infected by the bacteria Burkholderia at prescribed rates.
% Noh et al's "Facultative symbiont virulence determines horizontal .
% transmission rate without host specificity in Dictyostelium discoideum
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
rows = 20; %the width of the environment- each row increases width by 3.6 ...
                                                                        %cm
columns = 20; % the height of the environment- each row increases height...
                                                                %by 3.6 cm

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

%This creates the same "random" output everytime and the seed can be
%changed to create a new set of data
rng(34)

% set simulation's duration and variables:

%This is the length that the simulation will run for. Each iteration...
%represents an hour of time passing in the simulation
numIterations = 100;

%This is a parameter that can be used to turn integer visualization on and
%off. A value of 0 will visualize environments with only colored squares. A
%value of 1 will visualize environments with overlayed integer values for
%infected number and cluster size- for better readability. 
integerVisualization=1;

neighborhoodSize = 8; % size of neighborhood. In the cellular automata...
%model, this is a moore neighborhoods mean that the 8 adjacent cells will 
%be registered when updating the current cell. 

%this sets the interval at which amoebas will reproduce. Standard
%biological value would indicate that amoebas reproduce every 4 hours, so
%that value is set here. This can be changed depending on environmental
%conditions. 
reproductionTime=4; 

%This value determines the interval at which infection will spread in the
%population. Given that Burkholderia can reproduce around every 20 minutes,
%we anticipate that infection would spread horizontally  in adjacent 
%(clustered) amoebas on every iteration of the simulation. (Every 1 hour)
infectionTime=4;%every 4 hours infection will spread in the clusters

%At a certain point, infected amoeba will succumb to their infection and
%die. This value indicates how frequently in the simulation that happens.
%Can be changed by domain experts to increase biological accuracy. 
infectionLimit=12; %after 12 hours of persistent infection, amoeba will die

%There is biological evidence per Professor Suegene Noh that different
%reproduction rates may exist for infected and uninfected amoeba. This is
%implemented here. These values can be changed by domain experts to
%increase biological accuracy. 
infectedReproductionRate = 1.4;%how fast amoebas reproduce(1.2 = 20% growth,
%1.0 = 0% growth).
uninfectedReproductionRate=2;   

%This indicates the percentage of amoebas that will die off when we reach
%the infection limit. 25-35% of amoebas are believed to persist with
%infection so we assume that around 65% of them will succumb to their
%infection and die. This value can be changed by the domain expert to
%increase biological accuracy. 
infectionDeathRate=.65;

%this is the starting number of amoebas in the simulation.  
initialNumAmoebas = 1000; 

%This is the number of Amoebas in the simulation needed to create a
%fruiting body.
slugTotal = 300;

%Given the magnitude of amoebas in a simulation, agent based modeling for
%each individual amoeba would be computationally intensive. As such,
%aggregate modeling was chosen in order to better visualize large numbers
%of these unicellular organism. 
amoebasPerCell = 10; %number of amoebes per grid square (how we want to 
%initially divide up our amoeba population into aggregate "clusters"

%number of clusters started with
numClusters = initialNumAmoebas/amoebasPerCell;  

%sets a storage location for cluster positions
clusterPositions = [];  
                                            
%sets what proportion of amoebas are infected 0.0-1.0
percentInfected = 0.2; 

% Calculate how many infeced should be in each cell
infectedPerCell = round(amoebasPerCell * percentInfected);                                                      
                                                      
numInfected=numClusters*infectedPerCell; %use the initial clusters  
%infected to calculated the number of infected amoeba at the beginning of  
%thesimulation. 

percentSpore = 0.05; %The percent of amoebas that should be included in a 
% single spore when the fruiting body releases them
 
environment =  repmat({zeros(1,2)},rows,columns); % natural environment
%this consists of a cell array containing cells with 2 values in them.
%These will later be updated to store infected amoebas and cluster sizes


extEnvironment = repmat({zeros(1,2)},rows+2,columns+2); % environment with
%bounds part of the cellular automata process to ensure edge cases can read
%in neighbors and update. 

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
food = 4000; % Starting amount of food in the environment
starvationThreshold = 150; % Food level at which clusters start to clump
foodDecayRate = 1;%rate at which amoebas consume available food
foodList = 1:numIterations; %stores food availability at each point in 
                                                        %the simulation
                                                        
foodList(1) = food; %sets the first stored food value to the initial 
%allotted food for the environment



%array to keep track of iterations where fruiting body and spore dispersal
%happens for visualization efforts.
fruitingBodyIter = zeros(1,numIterations);


% initialize position and size of amoeba clusters
for i = 1:numClusters
    
    %sets the clusters position randomly in the envrionment
    clusterPos = [randi([1 rows]) randi([1 columns])];
    
    %stores randomly generated cluster positions if not already an existing
    %cluster at that position
    if isempty(clusterPositions)
        clusterPositions = [clusterPos];
    elseif sum(ismember(clusterPositions,clusterPos,'rows'),'all') <1   
        clusterPositions = [clusterPositions; clusterPos];
    end
    
    %Adds amoebas to the current ammount of amoebas in the cell it is 
    %randomly generating into.
    size = getClusterSize(clusterPos(1), clusterPos(2),environment);
    infection = getInfected(clusterPos(1), clusterPos(2),environment);
    environment{clusterPos(1), clusterPos(2)} = ...
     [size + amoebasPerCell,infection + infectedPerCell];

end

%keeps track of the number of clusters that are active in the environment.
% All clusters are initially active.
aliveClusters = length(clusterPositions);  

%update extended environment to contain the newly initiated environment
extEnvironment(2:rows+1, 2:columns+1) = environment;

%set first index of each storage list correctly
environmentList{1} = environment;
extEnvironmentList{1} = extEnvironment;
clusterPosList{1} = clusterPositions; 

%This function reads in a position value (x,y) and adds together the
%values of the adjacent moore neighborhood grid locations that are accessed
%by indexing from (-1,-1) to (1,1) and all the combinations in between .
%In the environment, empty 3.6 squarecentimeter grids that do not contain 
%amoebas have a clusterSize value of 0.
sumNeighbors = @(x, y, extEnvironment) (extEnvironment{x+1-1, y+1-1}(1) + ...
    extEnvironment{x+1-1, y+1}(1) + extEnvironment{x+1-1, y+1+1}(1) + ...
    extEnvironment{x+1, y+1-1}(1) + extEnvironment{x+1, y+1+1}(1) + ...
    extEnvironment{x+1+1, y+1-1}(1) + extEnvironment{x+1+1, y+1}(1) + ...
    extEnvironment{x+1+1, y+1+1}(1));


%Indexs for moving Amoeba clusters. These indexes represent the 8 index
%moves made to access the 8 neighbors in a moore neighborhood. 
indexMapping ={[-1 -1], [0 -1], [1 -1], [-1 0], [1 0], [-1 1],[0 1],[1 1]};

%create total counters to keep track of total amount of amoebas for
%visualization through simulation
amoebaTotals= zeros(1,numIterations);
amoebaTotals(1) = initialNumAmoebas;

%holding all total infection for simulation with each 
infectionTotals = zeros(1,numIterations);
infectionTotals(1) = numInfected;

%arrays to keep track of horizontal and vertical transmission for
%visualization
verticalTotals= zeros(1,numIterations);
horizontalTotals = zeros(1,numIterations);
 
% This loop implements our model for amoeba social behavior and
% Burkholderia infection by updating environment states based on prescribed
% logic, and storing each environment in a list that will later be used to
% visualize the cellular automata based simulation. 
for i = 2:numIterations
    %set environments properly for current iteration
    environment = environmentList{i-1}; %indicate which environment we will
    %be using to update the simulation for the next iteration
    
    extEnvironment = extEnvironmentList{i-1};%do the same for the extended
    %environment
    
    %identify the most recent cluster position list and set a variable to
    %make a new one based on movements that occur during this step in the
    %simulation
    oldClusterPos = clusterPosList{i-1};
    newClusterPos = [];  
    
    %counters to keep track of horizontal and vertical infections across
    %each iteration
    newHorizontalInfec = 0;
    newVerticalInfec = 0;

    %cycle through each cluster and move them accordingly based on
    %environmental factors
    for j = 1:aliveClusters
        
        %set position of the current cluster for updating by accessing that
        %cluster from the cluster position list
        clusterPos = oldClusterPos(j,:);
 
        %get current characteristics of cluster
        clusterSize = getClusterSize(clusterPos(1),clusterPos(2),...
            environment);
        clusterInfected = getInfected(clusterPos(1),clusterPos(2),...
            environment);
        clusterUninfected = clusterSize-clusterInfected;
        
        %Check to see if a cluster is big enough to form a slug
        if (fruitingBodyIter(i-1) ~= 1)&&...
            (clusterSize>= slugTotal)
            %Reset environment and positions to remvove the remaining
            %clusters, as only one slug will form.
            environment = repmat({zeros(1,2)},rows,columns);
            newClusterPos = clusterPos;
            
            %put slug back into the environment
            environment = setClusterSize(clusterPos(1),clusterPos(2),...
            clusterSize,environment);
            environment = setInfected(clusterPos(1),clusterPos(2),...
            clusterInfected,environment);
        
            %reset how many clusters are alive
            aliveClusters = 1;
            numAmoebas = clusterSize;
            
            %used for visualization to know when to add caption and to
            %ensure that spor dispersal happens
            fruitingBodyIter(i) = 1;
             
        %check to see if fruiting body has occured and then it is time for
        %spore dispersal.
        elseif fruitingBodyIter(i-1) == 1
            %find how amoebas per spore there are base on given parameters
            %above. Multiplied by 0.8 becaus 20% of amoebas die in the 
            %stalk of the fruiting body
            %round all calculations to end up with full amoebas and spores
            clusterSize = 0.8*clusterSize;
            clusterInfected = 0.8*clusterInfected;
            amoebasPerCell = round(clusterSize * percentSpore);
            infectedPerCell = round(clusterInfected* percentSpore);
            numSpores = round(clusterSize/ amoebasPerCell);
            
            %reset environment to remove the fruiting body
            environment = repmat({zeros(1,2)},rows,columns);
            newClusterPos = [];
            
            %reset food level
            foodList(i-1) = food;
            
            %reassign how many clusters and Amoebas will be alive after
            aliveClusters = numSpores;
            numAmoebas = amoebasPerCell*numSpores;
                %dispers spores throughout the environment.
                for s = 1:numSpores
                %sets the spore position randomly in the envrionment
                clusterPos = [randi([1 rows]) randi([1 columns])];
                
                %check to see if any spores exist at the randomly genrated 
                %placement, if not add it to cluster list
                if isempty(newClusterPos)
                    newClusterPos = [clusterPos];
                elseif sum(ismember(newClusterPos,clusterPos,'rows'),'all') <1   
                    newClusterPos = [newClusterPos; clusterPos];
                else
                    %if it is combinging into another spore remove it
                    aliveClusters = aliveClusters-1;
                end

                %get spores that are in the current cell of the environment 
                size = getClusterSize(clusterPos(1), clusterPos(2),...
                    environment);
                infection = getInfected(clusterPos(1), clusterPos(2),...
                    environment);


                %set number of infected clusters, all amoebas in cluster 
                %will be infected
                environment = setClusterSize(clusterPos(1),clusterPos(2),...
                        size + amoebasPerCell,environment);
                environment = setInfected(clusterPos(1),clusterPos(2),...
                        infection + infectedPerCell,environment);
                end
                
                %set an indicator for 
                fruitingBodyIter(i) = 2;
        else
            %set grid element in environment to be zero to show amoeabas 
            %have moved, this removes amoebas from their previous location
            environment{clusterPos(1),clusterPos(2)}= [0,0];
            extEnvironment(2:rows+1, 2:columns+1) = environment;

            %reproduce on reproductionTime- if it has been the number of hours
            %defined by reproductionTime, then the clusters will reproduce.
            %Unicellular reproduction only happens prior to and after slug
            %formation happens- there is no reproduction once amoebas have
            %fully clustered into a slug. 
            if rem(i,reproductionTime)==0 && foodList(i-1)>starvationThreshold
                %check if it is a reproduction iteration

                %remove current cluster size from total amount of amoebas
                numAmoebas = numAmoebas - clusterSize;

                %this increases the number of infected amoebas based on the
                %prescribed infectedReproductionRate. 
                clusterInfected = round(infectedReproductionRate ...
                    *clusterInfected);

                %this increases the number of uninfected amoebas based on the
                %prescribed unInfectedReproductionRate. 
                clusterUninfected=round(uninfectedReproductionRate*...
                    clusterUninfected);   

                %This recalculates the cluster size by adding together the 2
                %reproduced popualtions
                clusterSize = clusterInfected+clusterUninfected;

                %add updated cluster size to the total amount of amoebas
                numAmoebas = numAmoebas + clusterSize;
                
                %update vertical infection
                newVerticalInfec = newVerticalInfec + clusterInfected-...
                    clusterInfected/infectedReproductionRate;

            end

            %transmit infection on infectiontime. If it has been the number of
            %hours prescribed by infectionTime (i.e. infection spreads every
            %hour) then infection within clusters will be transmitted
            %horizontally. This mimics how in horizontal infection, amoebas in
            %close proximity to one another can infect one another. In this
            %case, clusters represent a certain number of amoebas in a 3.6
            %square centimeter space in the environment. As such, it makes
            %sense that those occupying the same space in the environment would
            %be the ones to infect one another. 
            if rem(i,infectionTime)==0 && clusterSize~=0 %this checks to make
                %sure it is an iteration that infection transmission happens on

               %calculate the percentage of the cluster that is currently 
               %infected
               infectedPercentage=(clusterInfected/clusterSize); 

               %there is biological evidence to show that infection spreads
               %with a linear relationship to the amount of infection in
               %the surrounding area. I.e. if more of your surrounding
               %environmentis infected, you have a higher chance of getting 
               %infected.
               %Thus, we calculate the number of amoebas that get infected
               %proportional to the number of amoebas that are infected in 
               %the surrounding environment (cluster). The same percentage 
               %of the unininfected amoebas becomes infected as was orignially 
               %infected in the cluster. 
               
               %keep track of new infection
               newInfection = 0;
               
               %go through each uninfected amoeba and calcuate a random 
               % number between 0 and 1, if less than infected percentage, 
               %then that amoeba becomes infected
               for amoeba = 1:clusterUninfected
                    probabilityOfInfection=rand;
                    if probabilityOfInfection<infectedPercentage
                      newInfection = newInfection+1;
                    end
               end

               %update cluster infection and unInfected
               clusterUninfected = clusterUninfected - newInfection;
               clusterInfected=clusterInfected + newInfection;
               
               %update total amount of horizontal infection happening
               newHorizontalInfec = newHorizontalInfec + newInfection;

            end

            %25-35% of infected amoebas persist with infection so around .65%
            %will die when they reach the infection limit iteration
            %checks to see if it is the prescribed interval in the simulation
            %where infected amoebas reach their limit and die off
            if rem(i,infectionLimit)==0

              %calculates the number of amoebas that die from infection based
              %on the infection death rate, round to be integer
              infectionKilled=round(infectionDeathRate*clusterInfected);

              %updates the cluster infected size to reflect the infected
              %amoebas that have died
              clusterInfected=clusterInfected-infectionKilled;

              %updates the cluster size to reflect the infected amoebas that
              %have died
              clusterSize=clusterSize-infectionKilled;

              %remove amoebas killed from total amount of amoebas
              numAmoebas = numAmoebas - infectionKilled;

            end

            %if there are neighbors that are ameobas, and the environment is 
            %is below the starvation threshold then combine with those 
            %neighbors
            %Empty cells have a value of 0, if sumNeighbors>0 it means
            %amoebas occupy adjacent cells. If the food has dipped below 
            %thestarvation threshold, this means social behavior is 
            %occurring and the amoebas will cluster together
            if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) >...
                    0 && foodList(i-1) < starvationThreshold

                %get list of neighbors that is shaped into 1d array
                neighbors = getNeighborSizes(clusterPos(1),clusterPos(2),...
                    extEnvironment);
                
                % find which direction to go and size of largest neighbor
                [maxNeighborSize, index] = max(neighbors);
                clusterMove = indexMapping{index};

                %get size of the neighboring cluster
                neighborPos= [clusterPos(1)+clusterMove(1),clusterPos(2)+...
                    clusterMove(2)];

                %update the clusterSize and number of infected once combined-
                %this means that the new cluster will combine the sizes and
                %infected amounts from the 2 clusters before they combined
                clusterSize = maxNeighborSize + clusterSize;
                clusterInfected =clusterInfected+getInfected(neighborPos(1),...
                    neighborPos(2),environment);

                %update the neighboring cell location in the environment to 
                %show the combined size and infection
                 environment = setClusterSize(neighborPos(1),...
                    neighborPos(2),clusterSize,environment);
                environment = setInfected(neighborPos(1),...
                    neighborPos(2),clusterInfected,environment);
                %update the external environment to contain new environment
                extEnvironment(2:rows+1, 2:columns+1) = environment;

                %since two of the clusters have combined, we have one less
                %active and updating cluster in the environment
                aliveClusters = aliveClusters - 1;
            else

                if foodList(i-1) < starvationThreshold && aliveClusters > 1
                    % Move towards the closest cluster in the environment,
                    % amoebas are attempting to group together into a slug due
                    % to their starving condition
                    clusterMove = findClosestCluster(clusterPos, oldClusterPos);
                else
                    %randomly choose movement for the cluster- amoebas continue
                    %to randomly move around the environment if there is no
                    %starvation and they are not clustering
                    indsa = randi([1, 8]);
                    clusterMove = indexMapping{indsa};
                end

                %reset the cluster position to reflect the movement choice made
                %above
                clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows),...
                    minRows), max(min(clusterPos(2) + clusterMove(2), columns), ...
                    minCols)];   
                
                %update new position in position list, if it isn't already
                %in the list of positions
                if isempty(newClusterPos)
                    newClusterPos = [clusterPos];
                elseif sum(ismember(newClusterPos,clusterPos,'rows'),'all') <1   
                    newClusterPos = [newClusterPos; clusterPos];
                else
                    %if it is combinging into another cluster remove it
                    aliveClusters = aliveClusters-1;
                end


                %update infection and size in the environment based on any 
                %updates made above due to the time in the simulation
                %Get any amoebas that may currently exist in the cell
                %that is being moved into
                size = getClusterSize(clusterPos(1), clusterPos(2),...
                    environment);
                infection = getInfected(clusterPos(1), clusterPos(2),...
                    environment);
                
                %update environment with correct size and infection
                environment = setClusterSize(clusterPos(1),clusterPos(2),...
                    clusterSize+size,environment);
                environment = setInfected(clusterPos(1),clusterPos(2),...
                    clusterInfected+infection,environment);
                %update the external environment to contain new environment
                extEnvironment(2:rows+1, 2:columns+1) = environment;
            end
        end
    end
        
    %get total number of infection and place it into array
    totalInfection = sum(getInfectedEnvironment(environment),"all");
    infectionTotals(i) = totalInfection;
    
    %get total number of amoebas and place it into array
    numAmoebas = sum(getSizeEnvironment(environment),"all");
    amoebaTotals(i) = numAmoebas;
    
    %put total amounts of horizontal and vertical infection into arrays
    verticalTotals(i)= newVerticalInfec;
    horizontalTotals(i) = newHorizontalInfec;
    
    %update current food and food variable list
    %food decreases based on the number of amoebas in the environment
    %and ensure that food doesn't go below zero. 
    updateFood = round(foodList(i-1) - (numAmoebas * foodDecayRate));
    foodList(i) = max([updateFood 0]);
    
    
    %update environments lists
    environmentList{i} = environment;
    extEnvironmentList{i} = extEnvironment;
    clusterPosList{i} = newClusterPos;

end
 
%Call  function to create graphs to see infection totals and v.s. total
%amoebas in the environment
graphInfectionTotals(infectionTotals,amoebaTotals, infectedReproductionRate,...
    fruitingBodyIter)
%Call  function to create graphs to see horizontal infection totals v.s. 
%vertical infection totals.
graphInfectionType(horizontalTotals, verticalTotals,...
    infectedReproductionRate)
%call the function show_CA_list to visualize the simulation based on stored
%values of all of the environments that were generated along with
%environment paramters, food information, whether or not integers should be
%visualized, the number of amoebas that we set and the number of infected
%amoebas
show_CA_List(environmentList,fruitingBodyIter, ...
      rows,columns,1, foodList,integerVisualization,initialNumAmoebas,numInfected);


%getter methods for cluster size and infection. returns either
%the clusterSize or the clusterInfected based on the given position in the
%given environment. 
function clusterSize = getClusterSize(row,col, environment)
    clusterSize = environment{row,col}(1);
end

function clusterInfected = getInfected(row,col,environment)
    clusterInfected = environment{row,col}(2);
end

%this function gets the sizes of the neighbors in the moore neighborhood
%and returns a list of their sizes
function neighbors = getNeighborSizes(row,col,extEnvironment)
     matrixEnvironment = cell2mat(extEnvironment);
     neighbors = reshape(matrixEnvironment(row-1+1:...
                row+1+1,2*(col+1)-1-2:2:2*(col+1)-1+2),1,[]);
     neighbors(5) = [];
end


%this function sets a new cluster size for a given row, column, and new
%size in the environment
function updatedCellArray = setClusterSize(row,col, size, environment)
    environment{row,col}(1) = size;
    updatedCellArray = environment;
end


%this function sets the infected value in the environment based on a given
%location, infected value and environment to update
function updatedCellArray = setInfected(row,col,infected, environment)
    environment{row,col}(2) = infected;
    updatedCellArray = environment;
end


%given that all of the values are stored in the same environment in a cell
%array, this function returns only the values of the infection in the
%environment in a matrix by turning the cell array into a matrix and only
%selecting the rows from the matrix that reflect infected values
function infectedEnvironment = getInfectedEnvironment(environment)
    matrixEnvironment = cell2mat(environment);
    [numRows,numCols] = size(matrixEnvironment);
    infectedEnvironment = matrixEnvironment(:,2:2:numCols);
end

%given that all of the values are stored in the same environment in a cell
%array, this function returns only the values of the cluster sizes in the
%environment in a matrix by turning the cell array into a matrix and only
%selecting the rows from the matrix that reflect size values
function sizeEnvironment = getSizeEnvironment(environment)
    matrixEnvironment = cell2mat(environment);
    [numRows,numCols] = size(matrixEnvironment);
    sizeEnvironment = matrixEnvironment(:,1:2:numCols-1);
end

%function to find closet cluster for starvation clustering. 
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

%This function is used to create a line plot of the total infection and
%total amoebas throughout the simulation.
function graphInfectionTotals(infectionTotals,amoebaTotals,...
    infectedReproductionRate,fruitingBodyIter)
    
    %create x axis
    time = 1:1:length(infectionTotals);
    %find where clustering happened
    verticalLinePositions = find(fruitingBodyIter);
    
    %plot new figure with lines
    figure(2)
    plot(time, amoebaTotals, 'b',time, infectionTotals, 'g');
    hold on;
    %plot all lines where clustering occured 
    for q = 1:2:length(verticalLinePositions)
        xpos = cast(verticalLinePositions(q),'uint8');
        plot([xpos xpos],[0 max(infectionTotals)],'-r','LineWidth', 0.8);
    end
    hold off;
    
    %add labels to graph
    xlabel('Time in Hours');
    ylabel('Total Number of Amoebas');
    legend('Total Number of Amoebas','Number of Amoebas Infected', ...
        'Fruiting Body Formation');
    title(['Total Number of Amoebas Over Time with Infected Reproduction Rate of ', ...
        num2str(infectedReproductionRate)]);
end

%This function is used to create two bar plots of the ammount of vertical 
%and horizontal infection that happens at each iteration.
function graphInfectionType(horizontalTotals,verticalTotals, ...
    infectedReproductionRate)
    
    %plot vertical infection over intervals
    time = 1:1:length(verticalTotals);
    figure(3)
    bar(time, horizontalTotals);
    xlabel('Time in Hours');
    ylabel('Total Number of Infected Amoebas');
    title(['Number of Horizontally Infected Amoebas Over Time with Infected Reproduction Rate of ', ...
         num2str(infectedReproductionRate)]);
    
    %plot horizontal infection over intervals
    figure(4)
    bar(time, verticalTotals);
    xlabel('Time in Hours');
    ylabel('Total Number of Infected Amoebas');
     title(['Number of Vertically Infected Amoebas Over Time with Infected Reproduction Rate of ', ...
         num2str(infectedReproductionRate)]);
end

%this function is used to visualize the simulation after environment values
%have been stored in the environmentList. It takes in the environmentList,
%environment parameters, the food list, whether or not to visualize
%integers, the number of amoebas at the start and the number of infected at
%the start. 
function [ ] = show_CA_List(environmentList, fruitingBodyIter,...
    rows,columns,interval, foodList,integerVisualization,numAmoebas,numInfected)
    figure(1)
    %this is used for colorbar visualization. determines the interval that
    %we will increment color visualization on by dividing the maximum value
    %by 20 and splitting it into sections
    infectedInterval=round(numInfected/20);
    sizeInterval=round(numAmoebas/20);
    disp(sizeInterval)
    
    %beggins iterating through each environment in the environment list
    
    for i=1:interval:length(environmentList)
        %set the environment to visualize
        environment = environmentList{i};
        
        %get the infected environment and the size environment using
        %functions described above
        infectedEnvironment=getInfectedEnvironment(environment);
        sizeEnvironment=getSizeEnvironment(environment);
        hold on;
        
        %we will visualize infection and clustering side by side in order
        %to better validate data and understand how infection spreads at
        %the same time that clustering social behavior occurs. This sets up
        %the subplot for infection
        subplot(1,2,1); 
        colormap(subplot(1,2,1),flipud(bone)); 
        imagesc(infectedEnvironment);
        caxis([0,numInfected]);
        colorbar;
        bar1=colorbar;
        
        %this increments the colorbar based on the scale of the initial
        %infection and the prescribed infectedInterval from above
        bar1.Ticks=[infectedInterval,2*infectedInterval,3*infectedInterval...
            ,4*infectedInterval,5*infectedInterval,6*infectedInterval,...
            7*infectedInterval,8*infectedInterval,9*infectedInterval,...
            10*infectedInterval,11*infectedInterval,12*infectedInterval...
            ,13*infectedInterval,14*infectedInterval,15*infectedInterval,...
            16*infectedInterval,17*infectedInterval,18*infectedInterval,...
            19*infectedInterval,20*infectedInterval,21*infectedInterval]; 
        
        %labels the colorbar based on the scale established above. 
        bar1.TickLabels={'empty',sprintf('%d infected amoeba',...
            infectedInterval),...
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
        
           %iterate through each cell in the infected environment
           for m = 1:size(infectedEnvironment)[1];
                for j= 1:size(infectedEnvironment)[2];

                    %if there is no infection in a given point in the
                    %environment, but there IS a cluster, visualize this
                    %cluster by putting an empty white rectangle on top of
                    %it
                    if (infectedEnvironment(m,j)==0 && sizeEnvironment(m,j)>0)
                        subplot(1,2,1);
                        rectangle('Position', [j-0.5,m-0.5, 1, 1], 'FaceColor',...
                        [1,1,1],'EdgeColor', 'k');
                        
                        %indicate that the number of infected amoebas in
                        %this cluster is 0
                        text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'red');
            
                    end
                end       
           end  
           
        %iterate through each cell in the infected environment
        for m = 1:size(infectedEnvironment, 1)
            for j = 1:size(infectedEnvironment, 2)
                
                    %if integer visualization is set to on, visualize the
                    %integer over each infected cluster that reflects the
                    %number of amoebas in the cluster that are infected
                    if (integerVisualization==1 && infectedEnvironment(m,j)~=0)
                    text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'red');
                    end
            end
        end
        
        hold off;
        
        %visualize figure and title it based on the frame in the
        %simulation, and the food available
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Infected Environment - Frame: %d, Food Available: %d'...
            , i, foodList(i)));
        set(gca,'YDir','reverse');

        
        hold on;
        
        %visualize cluster sizes as the second subplot
        colormap(subplot(1,2,2),parula);
        subplot(1,2,2);  
        imagesc(sizeEnvironment);
        caxis([0,numAmoebas]);
        colorbar;
        bar2=colorbar;
        
         %this increments the colorbar based on the scale of the initial
        %number of amoeas and the prescribed sizeInterval from above
        bar2.Ticks=[sizeInterval,2*sizeInterval,3*sizeInterval...
            ,4*sizeInterval,5*sizeInterval,6*sizeInterval,...
            7*sizeInterval,8*sizeInterval,9*sizeInterval,...
            10*sizeInterval,11*sizeInterval,12*sizeInterval...
            ,13*sizeInterval,14*sizeInterval,15*sizeInterval,...
            16*sizeInterval,17*sizeInterval,18*sizeInterval,...
            19*sizeInterval,20*sizeInterval,21*sizeInterval]; 
        
        %labels the colorbar based on the scale established above. 
        bar2.TickLabels={'empty',sprintf('%d infected amoeba',infectedInterval),...
            sprintf('%d total amoeba',2*sizeInterval),...
            sprintf('%d total amoeba',3*sizeInterval),...
            sprintf('%d total amoeba',4*sizeInterval),...
            sprintf('%d total amoeba',5*sizeInterval),...
            sprintf('%d total amoeba',6*sizeInterval),...
            sprintf('%d total amoeba',7*sizeInterval),...
            sprintf('%d total amoeba',8*sizeInterval),...
            sprintf('%d total amoeba',9*sizeInterval),...
            sprintf('%d total amoeba',10*sizeInterval),...
            sprintf('%d total amoeba',11*sizeInterval),...
            sprintf('%d total amoeba',12*sizeInterval),...
            sprintf('%d total amoeba',13*sizeInterval),...
            sprintf('%d total amoeba',14*sizeInterval),...
            sprintf('%d total amoeba',15*sizeInterval),...
            sprintf('%d total amoeba',16*sizeInterval),...
            sprintf('%d total amoeba',17*sizeInterval),...
            sprintf('%d total amoeba',18*sizeInterval),...
            sprintf('%d total amoeba',19*sizeInterval),...
            sprintf('%d + total amoeba',20*sizeInterval)};
        
       %if integer visualization is on, visualize an integer in text over
       %each cluster in the simulation to make it more readable and
       %understandable to viewers- shows how many amoebas are in each
       %visualized cluster
       for m = 1:size(sizeEnvironment, 1)
            for j = 1:size(sizeEnvironment, 2)
                if (integerVisualization==1 && sizeEnvironment(m,j)~=0)
                    text(j, m, num2str(sizeEnvironment(m, j)), ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'magenta');
                end
            end
        end
        
        
        
        hold off;
        
        %visualize figure
        axis equal; axis tight; axis xy;
        xlim([0.5, rows + 0.5]);
        ylim([0.5, columns + 0.5]);
        title(sprintf('Size Environment - Frame: %d, Food Available: %d', i,...
            foodList(i)));
        set(gca,'YDir','reverse');
        
        %This adds annotations if a fruiting body occurs or a spore dispers
        if fruitingBodyIter(i) == 1
            fba = annotation('textbox',[.2 .07 .75 .12], 'String',...
            ['This is a Fruiting Body of Amoebas. ',...
            'All other clusters removed from environment'],...
            'FitBoxToText','on', 'FontSize', 24);
        %This adds annotations if spore dispersal
        elseif fruitingBodyIter(i) == 2
            %remove previous annotation
            delete(fba);
            %add new annotation
            fba = annotation('textbox',[.2 .07 .75 .12], 'String',...
            'Spore dispersal, randomly placing spores in environment',...
            'FitBoxToText','on', 'FontSize', 24);
        elseif i>1 && fruitingBodyIter(i-1) == 2
            %remove previous annotation
            delete(fba);
        end

        %allow user to iterate through the visualized simulation by
        %pressing buttons
        w = waitforbuttonpress;
        
    end
end