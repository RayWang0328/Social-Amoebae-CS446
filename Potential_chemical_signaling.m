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

integerVisualization=1;
 
neighborhoodSize = 8; % size of neighborhood

reproductionTime=10; % every 10 frames amoebas will reproduce
reproductionRate = 2; %how fast amoebas reproduce(1.2 = 20% growth, 1.0 =
                        % 0% growth).

   
%MAY NEED SOME OF THESE DIFFUSION CONSTANTS FOR CHEMICAL SIGNALING
% r = 0.05; % Diffusion Constant
% coolingRate = 0.2; % Cooling Constant
numAmoebas = 1000; %number of total amoebas
amoebasPerCell = 10; %number of amoebes per section
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
clusterPositions = zeros(numClusters,2);

percentInfected = 0.2; %sets what proportion of amoebas are infected 0.0-1
numInfectedClusters = round(percentInfected*numClusters); %calculate number 
                                                      %of infected clusters
numInfected=percentInfected*numAmoebas;
 
environment =  repmat({zeros(1,3)},rows,columns); % natural environment
extEnvironment = repmat({zeros(1,3)},rows+2,columns+2); % environment with bounds
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



        diffusion_rate = 0.1;

        % Update chemical concentrations based on diffusion
        for r = 1:rows
            for c = 1:columns
             % Get chemical concentration in the current cell
                concentration = environment{r,c}(3);
                % Calculate new concentration based on diffusion from neighboring cells
                for dr = -1:1
                  for dc = -1:1
                      if r+dr >= 1 && r+dr <= rows && c+dc >= 1 && c+dc <= columns
                            % Calculate diffusion from neighboring cell
                            new_concentration = environment{r+dr,c+dc}(3) * diffusion_rate;
                            concentration = concentration + new_concentration;
                      end
                  end
                 end
        % Update the chemical concentration in the current cell
        environment{r,c}(3) = concentration;
            end
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
                environment{clusterpos(1), clusterpos(2)}(3) = environment{clusterpos(1), clusterpos(2)}(3) + 1;
                clusterMove = followchemicaltrail(clusterPos);
            else
                %randomly choose movement for the cluster
                indsa = randi([1, 8]);
                clusterMove = indexMapping{indsa};
            end
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows),...
                minRows), max(min(clusterPos(2) + clusterMove(2), columns), ...
                minCols)];          
                
            %update new position in position list
            
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
function closestMove = followchemicaltrail(clusterPos,environment)
    % Determine direction based on chemical concentration gradient
    max_concentration = -inf;
    best_direction = [0, 0];
    for dr = -1:1
        for dc = -1:1
            if dr == 0 && dc == 0
                continue; % Skip the current cell
            end
            neighbor_concentration = environment{clusterPos(1)+dr, clusterPos(2)+dc}(3);
            if neighbor_concentration > max_concentration
                max_concentration = neighbor_concentration;
                best_direction = [dr, dc];
            end
        end
    end
    % Move the amoeba in the direction of higher chemical concentration
    new_amoeba_row = amoeba_row + best_direction(1);
    new_amoeba_col = amoeba_col + best_direction(2);
    
    closestMove=[new_amoeba_row,new_amoeba_col];
end


function [ ] = show_CA_List(environmentList,...
    rows,columns,interval, foodList,integerVisualization,numAmoebas,numInfected)
    
    infectedInterval=numInfected/20;
    sizeInterval=numAmoebas/20
    
    for i=1:interval:length(environmentList)
        environment = environmentList{i};
        infectedEnvironment=getInfectedEnvironment(environment);
        sizeEnvironment=getSizeEnvironment(environment);
        hold on;
        
        
        subplot(1,2,1); 
        colormap(subplot(1,2,1),hsv); 
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
            sprintf('%d + infected amoeba',20*infectedInterval)};
        
           for m = 1:size(infectedEnvironment)[1]
                for j= 1:size(infectedEnvironment)[2]

            
                    if (infectedEnvironment(m,j)==0 && sizeEnvironment(m,j)>0)
                        subplot(1,2,1);
                        rectangle('Position', [j-0.5,m-0.5, 1, 1], 'FaceColor',...
                        [1,1,1],'EdgeColor', 'k');
                    
                        text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'k');
            
                    end
                end       
           end  
           
        for m = 1:size(infectedEnvironment, 1)
            for j = 1:size(infectedEnvironment, 2)
                    if (integerVisualization==1 && infectedEnvironment(m,j)~=0)
                    text(j, m, num2str(infectedEnvironment(m, j)),...
                        'HorizontalAlignment', 'center', 'VerticalAlignment'...
                        , 'middle', 'FontSize', 8, 'Color', 'k');
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