% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Fall, 2023
%

 
rows = 20; 
columns = 20;

% for extended grid checks to keep bounds
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 200;
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood

%MAY NEED SOME OF THESE DIFFUSION CONSTANTS FOR CHEMICAL SIGNALING
% r = 0.05; % Diffusion Constant
% coolingRate = 0.2; % Cooling Constant
numAmoebas = 50; %number of total amoebas
amoebasPerCell = 1; %number of amoebes per section
numClusters = numAmoebas/amoebasPerCell; %number of clusters started with
 
environment =  zeros(rows, columns); % natural environment
extEnvironment = zeros(rows+2, columns+2); % environment with bounds

%create lists to update information through out simulation
environmentList = zeros(rows, columns, numIterations);
extEnvironmentList = zeros(rows+2, columns+2, numIterations);
clusterPosList = zeros(numClusters,2,numIterations);  
clusterSizeList= zeros(numClusters,numIterations);




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
        clusterSizeList(j,i)=clusterSize;
        environment(clusterPos(1), clusterPos(2)) = 0;
        
        
        %if there are neighbors that are ameobas combine with them
        if sumNeighbors(clusterPos(1),clusterPos(2),extEnvironment) > 0
            
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
            
            clusterPosList(j,:,i)= clusterPos;
            
        end
        
    end
    
    
    
%     disp(environment)
%     w=waitforbuttonpress;
    environmentList(:,:,i) = environment;
    extEnvironmentList(:,:,i) = extEnvironment;

end
 
show_CA_List(environmentList, numAmoebas, clusterPosList,clusterSizeList,rows,columns,1);



    
function [ ] = show_CA_List(environmentList,numAmoebas,amoebasPosList,clusterSizeList,rows,columns,interval)

    for i=1:interval:length(environmentList)
        environment = environmentList(:,:,i);
        hold on;

       
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
        
        caxis([0,12]);
        lifeCycleColors=colorbar;
        lifeCycleColors.Ticks=[1,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5];   
        lifeCycleColors.TickLabels={'empty','1 amoeba','2 amoeba cluster','3 amoeba cluster',...
        '4 amoeba cluster','5 amoeba cluster','6 amoeba cluster','7 amoeba cluster','8 amoeba cluster','9 amoeba cluster','10+ amoeba cluster'};
         hold;
        
        
       
        % Plot mammals positions on top of the heat map
        for m = 1:numAmoebas
            
            if (1-clusterSizeList(m,i)*.1)>0
                greencolor=1-clusterSizeList(m,i)*.1;
      
            else greencolor=0; 
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

