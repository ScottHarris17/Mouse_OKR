function [upInputResistances, downInputResistances] = GetInputResistances()
% This function gets the input resistances of cells from cellList using the
% ramp stimulus.
% Rinput = %(V2-V1)/(I2-I1)... in other words, the slope of the ramp when
% plotted in the I-V plane

ups = grabFromFilter('UP Cells');
downs = grabFromFilter('DOWN Cells');

%choose which protocols to look at
priorityOrder_1 = {'Ramp'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; 'CurrentClamp'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {''};
cantHave_1 = [];
cantHave_2 = {'Cesium', 'Extracellular'}; %avoid anything with these tags
rigMandates = []; 

upInputResistances = [];
unames = {};
%loop through cell list
for i = 1:size(ups, 1)
    struct_i = ups{i, 2};
    %ask whether the current cell has the requested analysis
    [rampLoc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    if pass
        continue
    end
    unames(end+1) = {struct_i.cellID};

    inputRs = zeros(1, numel(rampLoc.Analysis_Results.EpochNumbers));
    for j = 1:numel(inputRs)
        restingPotential = rampLoc.Analysis_Results.restingMembranePotentialByEpoch(j); %mV
        thresholdPotential = rampLoc.Analysis_Results.thresholdByEpoch(j); %mV
        
        spikeTimes_j = rampLoc.Analysis_Results.allSpikeTimes{j};
        firstSpikeTime = spikeTimes_j(1) - 2; %ms... subtracting 2ms because spiketime indicates the peak of the spike
        preTime = rampLoc.meta.preTime; %ms
        
        timeToThreshold = firstSpikeTime - preTime; %ms - amount of time after stimulus onset that it took to reach threshold potential
        
        rampAmplitude = rampLoc.meta.rampAmplitude; %pA
        rampDuration = rampLoc.meta.stimTime; %ms
        
        currentSlope = rampAmplitude/rampDuration; %pA/ms
        
        currentAtThreshold = currentSlope * timeToThreshold; %pA - the amount of current that was being injected at threshold
        
        inputRs(j) = 1e9*(thresholdPotential-restingPotential)/currentAtThreshold; %Ohms
    end
    
    upInputResistances(end+1) = mean(inputRs);
end

        

downInputResistances = [];
dnames = {};
%loop through cell list
for i = 1:size(downs, 1)
    struct_i = downs{i, 2};
    %ask whether the current cell has the requested analysis
    [rampLoc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    if pass
        continue
    end
    dnames(end+1) = {struct_i.cellID};
    inputRs = zeros(1, numel(rampLoc.Analysis_Results.EpochNumbers));
    for j = 1:numel(inputRs)
        restingPotential = rampLoc.Analysis_Results.restingMembranePotentialByEpoch(j); %mV
        thresholdPotential = rampLoc.Analysis_Results.thresholdByEpoch(j); %mV
        
        spikeTimes_j = rampLoc.Analysis_Results.allSpikeTimes{j};
        firstSpikeTime = spikeTimes_j(1) - 2; %ms... subtracting 2ms because spiketime indicates the peak of the spike
        preTime = rampLoc.meta.preTime; %ms
        
        timeToThreshold = firstSpikeTime - preTime; %ms - amount of time after stimulus onset that it took to reach threshold potential
        
        rampAmplitude = rampLoc.meta.rampAmplitude; %pA
        rampDuration = rampLoc.meta.stimTime; %ms
        
        currentSlope = rampAmplitude/rampDuration; %pA/ms
        
        currentAtThreshold = currentSlope * timeToThreshold; %pA - the amount of current that was being injected at threshold
        
        inputRs(j) = 1e-3*(thresholdPotential-restingPotential)/(currentAtThreshold*1e-12); %Ohms
    end
    
    downInputResistances(end+1) = mean(inputRs);
end

end