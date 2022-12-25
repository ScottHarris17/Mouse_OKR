ups = grabFromFilter('UP - Spikes');
downs = grabFromFilter('DOWN Cells');

windowTime = 40; %ms
stepSize = 5; %ms

%choose which protocols to look at
priorityOrder_1 = {'Flash'}; %first field name should be one of these. Priority is in given order
priorityOrder_1b = {'Moving_Bar'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; '_all'; 'Extracellular'}; %child field name should be one of these. Priority is in given order
priorityOrder_2b = {'_SP10'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'Extracellular'};
mustHave_2b = {'Extracellular', '_ff'};
cantHave_1 = [];
cantHave_2 = {'Cesium', 'Intracellular', 'Potassium'}; %avoid anything with these tags
rigMandates = []; %use only cells with LED spot diameter = 500

upMaxFlash = [];
upDSIs = [];
upTCAreas = [];
upNormedAreas = [];
for i = 1:size(ups, 1)
    struct_i = ups{i, 2};
    
    [loc_flash, pass_flash] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    [loc_bar, pass_bar] = getLoc(struct_i, priorityOrder_1b, mustHave_1, cantHave_1, priorityOrder_2b, mustHave_2b, cantHave_2, rigMandates);
    if pass_flash || pass_bar
        continue
    end
    
    %calculate max FR for bar:
    preTime_i = loc_flash.meta.preTime;
    stimTime_i = loc_flash.meta.stimTime;
    tailTime_i = loc_flash.meta.tailTime;
    
    totalTime = preTime_i + stimTime_i + tailTime_i;
    
    if totalTime ~= 3000
        continue
    end
    spikeTimes_i = loc_flash.Analysis_Results.allSpikeTimes;
    
    numBins = round(3000/stepSize);
    results_i = zeros(numel(spikeTimes_i), numBins);
    for j = 1:numel(spikeTimes_i)
        s = spikeTimes_i{j};
        firingRates = spikeTimesToFiringRates(s, 0, 3000, windowTime, stepSize);
        results_i(j, :) = firingRates * 1000;
    end
    
    if j == 1
        meanResults = results_i;
    else
        meanResults = mean(results_i);
    end
    upMaxFlash(end + 1) = max(meanResults);
    
    %Calculate metrics for bar
    barMetrics = SpikeMetricHelper(loc_bar);
    upDSIs(end + 1) = barMetrics.DSI;
    upTCAreas(end+1) = barMetrics.TuningCurveArea;
    upNormedAreas(end+1) = barMetrics.NormalizedTuningCurveArea;
end

downMaxFlash = [];
downDSIs = [];
downTCAreas = [];
downNormedAreas = [];
for i = 1:size(downs, 1)
    struct_i = downs{i, 2};
    
    [loc_flash, pass_flash] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    [loc_bar, pass_bar] = getLoc(struct_i, priorityOrder_1b, mustHave_1, cantHave_1, priorityOrder_2b, mustHave_2b, cantHave_2, rigMandates);
    if pass_flash || pass_bar
        continue
    end
    
    %calculate max FR for bar:
    preTime_i = loc_flash.meta.preTime;
    stimTime_i = loc_flash.meta.stimTime;
    tailTime_i = loc_flash.meta.tailTime;
    
    totalTime = preTime_i + stimTime_i + tailTime_i;
    
    if totalTime ~= 3000
        continue
    end
    spikeTimes_i = loc_flash.Analysis_Results.allSpikeTimes;
    
    numBins = round(3000/stepSize);
    results_i = zeros(numel(spikeTimes_i), numBins);
    for j = 1:numel(spikeTimes_i)
        s = spikeTimes_i{j};
        firingRates = spikeTimesToFiringRates(s, 0, 3000, windowTime, stepSize);
        results_i(j, :) = firingRates * 1000;
    end
    
    if j == 1
        meanResults = results_i;
    else
        meanResults = mean(results_i);
    end
    downMaxFlash(end + 1) = max(meanResults);
    
    %Calculate metrics for bar
    barMetrics = SpikeMetricHelper(loc_bar);
    downDSIs(end + 1) = barMetrics.DSI;
    downTCAreas(end+1) = barMetrics.TuningCurveArea;
    downNormedAreas(end+1) = barMetrics.NormalizedTuningCurveArea;
end


%%
disp('TC Area')
[r_up, p_up] = corr(upMaxFlash', upTCAreas', 'type', 'spearman')
[r_down, p_down] = corr(downMaxFlash', downTCAreas', 'type', 'spearman')
disp('DSI:')
[r_up, p_up] = corr(upMaxFlash', upDSIs', 'type', 'spearman')
[r_down, p_down] = corr(downMaxFlash', downDSIs', 'type', 'spearman')
disp('Normalized Area')
[r_up, p_up] = corr(upMaxFlash', upNormedAreas', 'type', 'spearman')
[r_down, p_down] = corr(downMaxFlash', downNormedAreas', 'type', 'spearman')
