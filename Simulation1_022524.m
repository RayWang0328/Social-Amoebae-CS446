% Emmett Smith, Ray Wang, MJ Pennington
% CS346 
% Fall, 2023
%

 
rows = 20; 
columns = 20;
% for extended grid checks
minRows = 2;
minCols = 2;

% set simulation's duration and variables
numIterations = 100;
x = 1:rows;
y = 1:columns;
 
neighborhoodSize = 8; % size of neighborhood

%MAY NEED SOME OF THESE DIFFUSION CONSTANTS
r = 0.05; % Diffusion Constant
coolingRate = 0.2; % Cooling Constant
numAmoebas = 10; %number of mammals

 
environment =  zeros(rows, columns); % natural environment
extEnviorment = zeros(rows+2, columns+2); % environment with bounds

%create lists to update information through out simulation
environmentList = zeros(rows, columns, numIterations);
extEnviormentList = zeros(rows+2, columns+2, numIterations);
amoebasPosList = zeros(numAmoebas,2,numIterations); %rows 1..n-1 mammals,row n viper                                 


% initial positions of amoebas
for i = 1:numAmoebas
    amoebasPosList(i,:,1) = [randi([1 rows]) randi([1 columns])];
end


%set amoebas
for i = 1:numAmoebas
    amoebasPos = amoebasPosList(i,:,1);
    environment(amoebasPos(1), amoebasPos(2)) = 1;
end

 
 
 
% random movement function
randomMove = @(pos)[max(min(pos(1)+randi([-1,1]), rows), minRows) ...
    max(min(pos(2)+randi([-1,1]), columns), minCols)];
 
% viper movement towards mammal
moveTowards= @(targetPos, currentPos)(currentPos + sign(targetPos - currentPos));


%set first index of each list correctly
environmentList(:,:,1) = environment;
extEnviormentList(:,:,1) = extEnviorment;
 
% simulation loop
for i = 2:numIterations
    %set enviorments properly
    environment = environmentList(:,:,i-1);
    extEnviorment = extEnviormentList(:,:,i-1);
    
    %Indexs for moving Amoebas
    indexMapping ={[-1 -1], [0 -1], [1 -1], [-1 0], [1 0], [-1 1],...
            [0 1],[1 1]};

    
    for j = 1:numAmoebas
        amoebaPos = amoebasPosList(j,:,i-1);
        disp(amoebaPos)
        extEnviorment(amoebaPos(1), amoebaPos(2)) = 0;
        indsa = randi([1, 8]);
        amoebaMove = indexMapping{indsa};
        amoebaPos = [max(min(amoebaPos(1) + amoebaMove(1), rows), minRows), ...
                    max(min(amoebaPos(2) + amoebaMove(2), columns), minCols)];
                
        amoebasPosList(j,:,i)= amoebaPos;
    end
    disp(amoebasPosList(j,:,i));
   
end
 
show_CA_List(environmentList, numAmoebas, amoebasPosList,rows,columns,1);

    
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

        pause(0.1); % wait for a moment to proceed to the next frame
        
    end
end

