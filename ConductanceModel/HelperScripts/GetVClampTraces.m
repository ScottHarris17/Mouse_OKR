function [excitationMeans, inhibitionMeans] = GetVClampTraces(useUpCells, useDownCells)
%Returns mean VClamp Current traces for each direction (excitation and inhibition).
%User can choose whether to simulate upCells, downCells, or both. Ex: set
%useUpCells input value to 1 to use data from up cells and 0 to not use data
%from up cells. upCells and downCells cannot both be zero

if useUpCells + useDownCells == 0
    errormsg('Must include data from at least one cell type')
    return
end

orientations = [0:45:315];

%choose which protocols to look at
priorityOrder_1 = {'Moving_Bar'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; 'Intracellular'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'Intracellular', 'Excitation'};
mustHave_2b = {'Intracellular', 'Inhibition'};
cantHave_1 = {'Var'};
cantHave_2 = {'Extracellular', 'CurrentClamp'}; %avoid anything with these tags
rigMandates = []; %use only cells with LED spot diameter = 500

excitationAll = cell(1, 8); %initialize holders
inhibitionAll = cell(1, 8); %initialize holders
%ups first

if useUpCells
    ups = grabFromFilter('UP Cells');
    for i = 1:size(ups, 1)
        struct_i = ups{i, 2};
        
        [loc_excitation, pass1] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
        [loc_inhibition, pass2] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2b, cantHave_2, rigMandates);
        
        if ~pass1
            loc = loc_excitation;
            for j = 1:numel(orientations)
                try
                index_j = loc.Analysis_Results.meanTracesByOrientation.orientationOrder == orientations(j);
                catch
                    me
                end
                trace = loc.Analysis_Results.meanTracesByOrientation.meanTraces(index_j, :);
                soFar = excitationAll{j};
                plusThis = [soFar;trace];
                excitationAll{j} = plusThis;
            end
        end
        
        if ~pass2
            loc = loc_inhibition;
            for j = 1:numel(orientations)
                index_j = loc.Analysis_Results.meanTracesByOrientation.orientationOrder == orientations(j);
                trace = loc.Analysis_Results.meanTracesByOrientation.meanTraces(index_j, :);
                soFar = inhibitionAll{j};
                plusThis = [soFar;trace];
                inhibitionAll{j} = plusThis;
            end
        end
    end
end
%Repeat for downs
if useDownCells
    downs = grabFromFilter('DOWN Cells');
    downOrientations = [180:45:315, 0:45:135]; %rotate down responses so that can combine with up responses
    for i = 1:size(downs, 1)
        struct_i = downs{i, 2};
        
        [loc_excitation, pass1] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
        [loc_inhibition, pass2] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2b, cantHave_2, rigMandates);
        
        
        if ~pass1
            loc = loc_excitation;
            for j = 1:numel(orientations)
                index_j = loc.Analysis_Results.meanTracesByOrientation.orientationOrder == orientations(j);
                trace = loc.Analysis_Results.meanTracesByOrientation.meanTraces(index_j, :);
                soFar = excitationAll{j};
                plusThis = [soFar;trace];
                excitationAll{j} = plusThis;
            end
        end
        
        if ~pass2
            loc = loc_inhibition;
            for j = 1:numel(downOrientations)
                index_j = loc.Analysis_Results.meanTracesByOrientation.orientationOrder == downOrientations(j);
                trace = loc.Analysis_Results.meanTracesByOrientation.meanTraces(index_j, :);
                soFar = inhibitionAll{j};
                plusThis = [soFar;trace];
                inhibitionAll{j} = plusThis;
            end
        end
    end
end

%Get the excitation traces. Remember, excitation is supposed to be uniform
%across stimulus directions. Space clamp reduces recorded excitation.
%Therefore, the best approximation of excitation at every time point is the
%maximum of excitation across all stimulus directions
excitationMeans = zeros(numel(orientations), 5500);
for i = 1:numel(excitationAll)
    trace_j = excitationAll{i};
    thisTrace = mean(trace_j) - mean(trace_j(1:1000));%baseline subtract
    %thisTrace(1:1750) = 0; %option to zero the baseline
    excitationMeans(i, :) = thisTrace;
end

%because excitation should be directionally symmetric, scale all
%excitatory responses up to the maximum response (space clamp causes
%excitatory traces to look smaller than they should be in certain
%directions)
excitationMeanPeaks = min(excitationMeans, [], 2);
overallMin = min(excitationMeanPeaks);
scaleFactors = overallMin./excitationMeanPeaks;
for i = 1:size(excitationMeans, 1)
    excitationMeans(i, :) = excitationMeans(i, :) * scaleFactors(i);
end
%then take the maximum across all directions for each point to control for
%possible width changes.
maxExcitationTrace = min(excitationMeans);
maxExcitationTrace = maxExcitationTrace - mean(maxExcitationTrace(1:1000));
excitationMeans = repmat(maxExcitationTrace, numel(orientations), 1);

%get inhibition traces
inhibitionMeans = zeros(numel(orientations), 5500);
for i = 1:numel(inhibitionAll)
    trace_j = inhibitionAll{i};
    thisTrace = mean(trace_j) - mean(trace_j(1:1000));%baseline subtract
    %thisTrace(1:1750) = 0; %option to zero the baseline
    inhibitionMeans(i, :) = thisTrace;
end
end