%distance between spike threshold and resting membrane potential

ups = grabFromFilter('UP Cells');
downs = grabFromFilter('DOWN Cells');


%choose which protocols to look at
priorityOrder_1 = {'Flash'}; %first field name should be one of these. Priority is in given order
priorityOrder_1b = {'Ramp'}; 
priorityOrder_2 = {'_ff'; 'SP10'; 'CurrentClamp'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {''};
cantHave_1 = [];
cantHave_2 = {'Cesium', 'Extracellular'}; %avoid anything with these tags
rigMandates = []; %use only cells with LED spot diameter = 500

load cutoffValues.mat

upRestingAndThresholds = [];
downRestingAndThresholds = [];

for i = 1:size(cutoffValues, 1)
    
    cname = cutoffValues{i, 1};
    isUp = 0;
    isDown = 0;
    if ismember(cname, ups(:, 1)) || ismember(cname, downs(:, 1))
        [isUp, upIndex] = ismember(cname, ups(:, 1));
        [isDown, downIndex] = ismember(cname, downs(:, 1));
        if isUp
            struct_i = ups{upIndex, 2};
        else
            struct_i = downs{downIndex, 2};
        end

        [flashLoc, pass1] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
        [rampLoc, pass2] = getLoc(struct_i, priorityOrder_1b, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
        
        if pass1 || pass2
            continue
        end
        
        restingVm = flashLoc.Analysis_Results.restingMembranePotentialByEpoch(1);
        threshold = rampLoc.Analysis_Results.averageThreshold;
        
        if isUp
            upRestingAndThresholds = [upRestingAndThresholds, [restingVm; threshold; restingVm-threshold]];
        else
            downRestingAndThresholds = [downRestingAndThresholds, [restingVm; threshold; restingVm-threshold]];
        end
    end
end

%upRestingAndThresholds(3,7:8) = NaN;
upRestingAndThresholds(upRestingAndThresholds(3, :) > 0) = 0;
downRestingAndThresholds(downRestingAndThresholds(3, :) > 0) = 0;
edgs = [-30:5:10];

figure
title('Difference between Resting Vm and Threshold')
hold on
h1 = histogram(upRestingAndThresholds(3, :), 'BinEdges', edgs, 'FaceColor', 'k');
h2 = histogram(downRestingAndThresholds(3, :), 'BinEdges', edgs, 'FaceColor', 'm');

h2.BinEdges = h1.BinEdges;
p = ranksum(upRestingAndThresholds(3, :), downRestingAndThresholds(3, :))