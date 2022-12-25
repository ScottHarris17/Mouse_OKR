%Comparing holding currents for flash stimulus of UP and DOWN cells

ups = grabFromFilter('UP Cells');
downs = grabFromFilter('DOWN Cells');


%choose which protocols to look at
priorityOrder_1 = {'Flash'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; 'Intracellular'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'Intracellular', 'Excitation'};
cantHave_1 = [];
cantHave_2 = {'Extracellular', 'CurrentClamp'}; %avoid anything with these tags
rigMandates = []; %

upHoldingCurrent = [];
downHoldingCurrent = [];
upCoords = zeros(0, 2);
downCoords = zeros(0, 2);
%ups first
for i = 1:size(ups, 1)
    struct_i = ups{i, 2};
    
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    if pass
    	continue
    end
    
    holdingCurrent = loc.Analysis_Results.holdingCurrentByEpoch(1);
    if holdingCurrent < -200 %probably a bad recording
        continue
    end
    upHoldingCurrent(end+1) = holdingCurrent;
    upCoords(end+1, :) = struct_i.coordinates.polar;
end

%downs second
for i = 1:size(downs, 1)
    struct_i = downs{i, 2};
    
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    if pass
    	continue
    end

    holdingCurrent = loc.Analysis_Results.holdingCurrentByEpoch(1);
    if holdingCurrent < -200
        continue
    end
    downHoldingCurrent(end+1) = holdingCurrent;
    downCoords(end+1, :) = struct_i.coordinates.polar;
end

edgs = [-300:30:100];
figure
title('Holding Current')
hold on
h1 = histogram(upHoldingCurrent, 'BinEdges', edgs, 'FaceColor', 'black', 'Normalization','probability');
h2 = histogram(downHoldingCurrent, 'BinEdges', edgs, 'FaceColor', 'magenta', 'FaceAlpha', 0.6, 'Normalization', 'Probability');
h1.BinWidth = h2.BinWidth;
xlabel('pA')
ylabel('fraction')
legend('Ups', 'Downs')
box off

p_holdingCurrent = ranksum(upHoldingCurrent, downHoldingCurrent)

[horizontalComponentUps, verticalComponentUps] = ...
    pol2cart(upCoords(:, 1), upCoords(:, 2));
[horizontalComponentDowns, verticalComponentDowns] = ...
    pol2cart(downCoords(:, 1), downCoords(:, 2));

figure
title('Holding Current by Vertical Position')
hold on
scatter(verticalComponentUps, upHoldingCurrent, 'k', 'filled')
scatter(verticalComponentDowns, downHoldingCurrent, 'm', 'filled')
xlabel('Vertical Coordinate')
ylabel('Holding Current')

figure
title('Holding Current by Horizontal Position')
hold on
scatter(horizontalComponentUps, upHoldingCurrent, 'k', 'filled')
scatter(horizontalComponentDowns, downHoldingCurrent, 'm', 'filled')
xlabel('Horizontal Coordinate')
ylabel('Holding Current')

%% Compute what the corresponding difference in resting Vm might be:
InputResistanceUp = 2.74e+08; %average input resistance of cells from current clamp
InputResistanceDown = 1.815e+08;
%difference in millivolts between resting Vms of Up and Down cells:
RestingVm_Difference = 1e-9*(InputResistanceUp*mean(upHoldingCurrent) - ...
    InputResistanceDown*mean(downHoldingCurrent))

    
