% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Fall, 2023
%

 
rows = 5; 
columns = 5;

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 20;
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood

%MAY NEED SOME OF THESE DIFFUSION CONSTANTS FOR CHEMICAL SIGNALING
% r = 0.05; % Diffusion Constant
% coolingRate = 0.2; % Cooling Constant
numAmoebas = 5; %number of total amoebas
amoebasPerCell = 1; %number of amoebes per section
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
 
environment =  zeros(rows, columns); % natural environment
extEnvironment = zeros(rows+2, columns+2); % environment with bounds

%create lists to update information through out simulation
environmentList = zeros(rows, columns, numIterations);
extEnvironmentList = zeros(rows+2, columns+2, numIterations);
clusterPosList = zeros(numClusters,2,numIterations);                                 


% initialize position and size of amoeba clusters
for i = 1:numClusters
    clusterPos = [randi([1 rows]) randi([1 columns])];
    environment(clusterPos(1), clusterPos(2)) = amoebasPerCell;
    clusterPosList(i,:,1) = clusterPos;
end

extEnvironment(2:rows+1, 2:columns+1) = environment;
 
% COULD BE USED TO MOVE AMOEBAS TOWARDS OTHERS
% moveTowards= @(targetPos, currentPos)(currentPos + sign(targetPos - ...
%     currentPos));


%set first index of each list correctly
environmentList(:,:,1) = environment;
extEnvironmentList(:,:,1) = extEnvironment;

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
    
    %disp(environment)
   
    %cycle through each cluster and move them accordingly accordingly
    for j = 1:numClusters
        %set position and cluster size
        clusterPos = clusterPosList(j,:,i-1);

        if(clusterPos == [0,0])
             continue;
        end
        clusterSize = environment(clusterPos(1), clusterPos(2));
        disp(clusterSize)
        environment(clusterPos(1), clusterPos(2)) = 0;
        
        
        %if there are neighbors that are ameobas combine with them
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0
            disp("combined")
            %get neighbors and shape into 1d array
            neighbors = reshape(extEnvironment(clusterPos(1)+1-1:...
                clusterPos(1)+1+1,clusterPos(2)-1+1:clusterPos(2)+1+1)...
                ,1,[]);
            neighbors(5) = []; % remove the clusters original position
            
            [maxNeighborSize, index] = max(neighbors);% find which direction to go
            clusterMove = indexMapping{index};
            
            
            %get size of the neighboring cluster and combine the two
            environment(clusterPos(1)+clusterMove(1),clusterPos(2)+...
                clusterMove(2)) = maxNeighborSize + clusterSize;
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            
            %move comined cluster off the screen and remove from simulation 
             clusterPosList(j,:,i)= [0,0];
        else
            %randomly choose movement for the cluster
            indsa = randi([1, 8]);
            clusterMove = indexMapping{indsa};
            clusterPos = [max(min(clusterPos(1) + clusterMove(1), rows), minRows), ...
                        max(min(clusterPos(2) + clusterMove(2), columns), minCols)];
                    
            %update relevant variables
            environment(clusterPos(1), clusterPos(2)) = clusterSize; 
            extEnvironment(2:rows+1, 2:columns+1) = environment;
            disp(environment(clusterPos(1), clusterPos(2)))
            clusterPosList(j,:,i)= clusterPos;
            
        end
        
    end
    
%     disp(environment)
%     w=waitforbuttonpress;
    environmentList(:,:,i) = environment;
    extEnvironmentList(:,:,i) = extEnvironment;

end
 
show_CA_List(environmentList, numAmoebas, clusterPosList,rows,columns,1);

    
function [ ] = show_CA_List(environmentList,numAmoebas,amoebasPosList,rows,columns,interval)

    for i=1:interval:length(environmentList)
        environment = environmentList(:,:,i);
        hold on;

        % Display the heat map
        imagesc(environment);
        colormap('hot'); % use a colormap suitable for heatmaps
        colorbar; % optional: display a colorbar indicating temperature values
        
       
        % Plot mammals positions on top of the heat map
        for m = 1:numAmoebas
            rectangle('Position', [amoebasPosList(m,2,i)-0.5,...
                    amoebasPosList(m,1,i)-0.5, 1, 1], 'FaceColor', 'g', ...
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

