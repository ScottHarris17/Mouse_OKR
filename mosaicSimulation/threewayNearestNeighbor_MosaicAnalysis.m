function [resultsStruct, nn3fig1, nn3fig2] = threewayNearestNeighbor_MosaicAnalysis(centers, pixelsPerMM)

thresholdDistances = 0:10:400; %um

thresholdDistances_pix = thresholdDistances * pixelsPerMM/1000;

violationsPerCellByThreshold = cell(numel(thresholdDistances), 2);
vs = zeros(numel(thresholdDistances), 1);
totalPairs = zeros(numel(thresholdDistances), 1);

%% count total pairs and total fraction of pairs that are members of a triangle
for n = 1:numel(thresholdDistances)
    threshold_n = thresholdDistances_pix(n);
    totalComparisons_n = 0;
    violationsTotal = 0;
    totalPairs_n = 0;
    
    %first find all cells that are within threshold of the reference cell
    for i = 1:size(centers, 1)-2
        center_i = centers(i, :);
        inPlay = zeros(0, 2);
        for j = i+1:size(centers, 1) %start at i+1 so as not to double count pairs
            center_j = centers(j, :);
            if cartDist(center_i, center_j) < threshold_n
                inPlay = [center_j; inPlay]; %in play holds all cells that are within threshold of the reference cell
            end
        end
        
        totalPairs_n = totalPairs_n + size(inPlay, 1);
        %next, out of the previously identified cells, find how many are
        %also within threshold of each other
        for k = 1:size(inPlay, 1)-1
            center_k = inPlay(k, :);
            for m = 1:size(inPlay, 1)
                if k == m
                    totalComparisons_n = totalComparisons_n + 1;
                    continue
                end
                center_m = inPlay(m, :);
                if cartDist(center_k, center_m) < threshold_n
                    violationsTotal = violationsTotal + 1;
                    break %break once the first violation is found because we are counting the fraction of pairs that are also part of triples and don't care whether a pair is a member of multiple triangles here
                end
            end
        end 
    end
    totalPairs(n) = totalPairs_n;
    vs(n) = violationsTotal;
end


%% Loop again to count for each cell (this counts total number redundantly)
for n = 1:numel(thresholdDistances)
    threshold_n = thresholdDistances_pix(n);
    violationsPerCell = zeros(size(centers, 1), 1);
    
    %first find all cells that are within threshold of the reference cell
    for i = 1:size(centers, 1)
        violationsCell_i = 0;
        center_i = centers(i, :);
        inPlay = zeros(0, 2);
        for j = 1:size(centers, 1)
            center_j = centers(j, :);
            if cartDist(center_i, center_j) < threshold_n && j ~= i
                inPlay = [center_j; inPlay]; %in play holds all cells that are within threshold of the reference cell
            end
        end
        
        %next, out of the previously identified cells, find how many are
        %also within threshold of each other
        for k = 1:size(inPlay, 1)-1
            center_k = inPlay(k, :);
            for m = k+1:size(inPlay, 1)
                center_m = inPlay(m, :);
                if cartDist(center_k, center_m) < threshold_n
                    violationsCell_i = violationsCell_i + 1;
                end
            end
        end 
        violationsPerCell(i) = violationsCell_i;
    end
    
    violationsPerCellByThreshold{n, 1} = thresholdDistances(n);
    violationsPerCellByThreshold{n, 2} = violationsPerCell;
end

%% compile results and make figures

resultsStruct.centers = centers;
resultsStruct.pixPerMM = pixelsPerMM;
resultsStruct.thresholdsUM = thresholdDistances;
resultsStruct.totalPairsByThreshold = totalPairs;
resultsStruct.fractionOfPairsInThreesome = vs./totalPairs;
resultsStruct.violationsPerCellByThreshold = violationsPerCellByThreshold;

%figure 1 shows all of the violations at each distance as a fraction of the
%total number of neighbors at that distance
normalizedViolations = vs./totalPairs;
normalizedViolations(isnan(normalizedViolations)) = 0;
nn3fig1 = figure;
title('Fraction of Pairs that are Triplets')
hold on
plot(thresholdDistances, normalizedViolations')
xlabel('Distance From Reference Cells (um)')
ylabel('Fraction of Pairs')
hold off

%figure 2 is a version of the nearest neighbor histogram for each cell (smallest
%threshold size where a threesome is found)
threewayNearestDistances = zeros(size(centers,1)-2, 1);
for i = 1:size(centers, 1) %for each cell
    for j = 1:size(violationsPerCellByThreshold, 1) %for each threshold
        threshold_j = violationsPerCellByThreshold{j, 2};
        if threshold_j(i) > 1 %value of cell_i at threshold_j
            threewayNearestDistances(i) = violationsPerCellByThreshold{j, 1};
            break
        end
    end
end

threewayNearestDistances(threewayNearestDistances == 0) = nan; %convert nearest neighbor distance of 0 to nan;
resultsStruct.threewayNearestDistances = threewayNearestDistances;

nn3fig2 = figure;
histogram(threewayNearestDistances, 'Normalization', 'probability')

hold on
title('Nearest Threeway Neighbor Distances')
xlabel('Distance um')
ylabel('Fraction')
xlim([0 max(thresholdDistances)]);

end


    