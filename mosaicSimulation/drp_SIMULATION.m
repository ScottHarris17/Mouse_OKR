function meanDRP = drp_SIMULATION(centers, numBins, retinaRadius)

%% Compute density recovery profile given center coordinates and mask
radiiSteps = logspace(1, log10(retinaRadius*2), numBins); %in pixels

results = zeros(size(centers, 1), numel(radiiSteps));
for i = 1:size(centers, 1)    
    %first calculate a distance matrix
    centers_i = centers(i, :);
    
    %step through the radii sizes       
    for j = 1:numel(radiiSteps)        
        radiiSteps_j = radiiSteps(j);
        %calculate available area inside the mask for the current radius
        circs = [0, 0, retinaRadius/1000; centers_i./1000, radiiSteps_j/1000];
        OverlapMat = area_intersect_circle_analytical(circs);
        
        availableArea = OverlapMat(1, 2); %mm^2
        
        
        %step through the remaining cells and count how many are in the
        %available area. Assuming all cells are in the mask.
        cellsCounted = -1; %start at -1 b/c each cell will also count itself
        for k = 1:size(centers, 1)
            d = cartDist(centers_i, centers(k, :));
            if d < radiiSteps_j
                cellsCounted = cellsCounted + 1;
            end
        end
        density = cellsCounted/availableArea; %cells/mm^2
        
        results(i, j) = density;
    end        
end

meanDRP = mean(results);


end
