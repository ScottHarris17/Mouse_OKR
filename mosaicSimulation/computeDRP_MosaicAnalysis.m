function [resultsStruct, drpfig] = computeDRP_MosaicAnalysis(img, centers, mask, pixelsPerMM)
%% Compute density recovery profile given center coordinates and mask
retinaArea = sum(mask(:));

bw = rgb2gray(img);

numBins = 100;
radiiSteps = logspace(0, log10(max(size(bw))), numBins); %in pixels

numberOfSteps = size(centers, 1) * numel(radiiSteps);
results = zeros(numberOfSteps, 7);
count = 1;
disp(['Total cells to analyze = ', num2str(size(centers, 1))])
w = waitbar(0, 'Computing...');
for i = 1:size(centers, 1)
    tic
    waitbar(i/size(centers,1), w, 'Computing...')
    
    %first calculate a distance matrix
    distances = zeros(size(bw));
    centers_i = centers(i, :);
    parfor j = 1:size(bw, 1)
        distances_k = zeros(1, size(bw,2)); %have to do this weird allocation for parfor but it speeds it up significantly
        for k = 1:size(bw, 2)
            distances_k(k) = pdist([centers_i; j, k], 'euclidean');
        end
        distances(j, :) = distances_k; 
    end
    
    %step through the area radii sizes       
    for j = 1:numel(radiiSteps)        
        %calculate available area inside the mask for the current radius
        availableAreaPix = sum(sum(mask.*(distances < radiiSteps(j))));
        availableArea = availableAreaPix/(pixelsPerMM^2); %convert to sq mm
        
        %step through the remaining cells and count how many are in the
        %available area. Assuming all cells are in the mask.
        cellsCounted = -1; %start at -1 b/c each cell will also count itself
        uniqueCellsCounted = 0;
        for k = 1:size(centers, 1)
            d = cartDist(centers_i, centers(k, :));
            if d < radiiSteps(j)
                cellsCounted = cellsCounted + 1; %this number double counts cells
                if k > i
                    uniqueCellsCounted = uniqueCellsCounted + 1; %doesn't double count cells, but then can't necessarily get an accurate density measurement
                end
            end
        end
        
        density = cellsCounted/availableArea;
        
        adjustedRadius = sqrt(availableArea/pi); %this is the radius of the circle with area = availableArea
        
        results(count, :) = [i, radiiSteps(j)/pixelsPerMM, availableArea, adjustedRadius,...
            uniqueCellsCounted, cellsCounted, density];
        
        count = count + 1;
    end        
    disp(['Cell Number ', num2str(i)])
    toc
end
close(w)

colNames = {'Cell Number', 'Radius Size', 'Available Area', 'Adjusted Radius', 'Unique Cells', 'Total Cells', 'Density'};
resultsTable = array2table(results, 'VariableNames', colNames);

resultsStruct = struct();
resultsStruct.data.centers = centers;
resultsStruct.data.mask = mask;
resultsStruct.data.drp = resultsTable;

%group the data into bins with equal radii
radii = results(:, 2);
[binNumbers, edges] = discretize(radii, radiiSteps/pixelsPerMM); %edges are in MM
densities = results(:, 7);

%keep track of the average and se density per bin
averagePerBin = zeros(1, numBins); %numBins is defined above. It is equal to length(radiiSteps)
SEPerBin = zeros(1, numBins);
for i = 1:numBins
    densityVals = [];
    for j = 1:numel(radii)
        if binNumbers(j) == i
            densityVals = [densityVals, densities(j)]; %have already computed densities above so now just tallying up the densities for each bin
        end
    end
    if numel(densityVals) == 0
        averagePerBin(i) = 0;
        SEPerBin(i) = 0;
        continue
    end
    averagePerBin(i) = mean(densityVals);
    SEPerBin(i) = std(densityVals)/sqrt(length(densityVals));
end
resultsStruct.pixPerMM = pixelsPerMM;
resultsStruct.Analysis.edges = edges;
resultsStruct.Analysis.binNumbers = binNumbers;
resultsStruct.Analysis.averagePerBin = averagePerBin;
resultsStruct.Analysis.SEPerBin = SEPerBin;

meanCellsPerMMsq = size(centers, 1)/(retinaArea/pixelsPerMM^2);

drpfig = figure;
title('Density Recovery Profile ')
hold on
box off
semilogx(edges*1000, averagePerBin)
set(gca, 'XScale', 'log');
xlabel('um')
ylabel('cells/mm^2')
plot(edges*1000, repmat(meanCellsPerMMsq, 1, numel(edges)), '--r');


end
