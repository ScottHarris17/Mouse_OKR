% This script compares the firing rates of UP and DOWN cells that came from
% the same retinas. Analyses here are based on firing rates in response to
% the flash stimulus.

downs = grabFromFilter('Downs');
ups = grabFromFilter('Ups');

%start by counting the number of UP and DOWN cells for each retina that
%you've recorded and store in 'retinas' matrix
retinas = [];
for i = 1:size(ups, 1)
    retina_i = ups{i, 1};
    name = str2double(retina_i(6:7));
    if ismember(name, retinas)
        retinas(retinas(:, 1) == name, 2) = retinas(retinas(:, 1) == name, 2) + 1;
    else
        retinas(end+1, :) = [name, 1, 0];
    end
end

for i = 1:size(downs, 1)
    retina_i = downs{i, 1};
    name = str2double(retina_i(6:7));
    if ismember(name, retinas)
        retinas(retinas(:, 1) == name, 3) = retinas(retinas(:, 1) == name, 3) + 1;
    else
        retinas(end+1, :) = [name, 0, 1];
    end
end

c_s = cell(0, 3);  %c_s will be a cell array that holds 3 things. Column 1 = animal name. Column 2 = up cell cell params struct. Column 3 = down cell cell param structs
for i = 1:size(retinas, 1)
    
    if retinas(i, 2) <2 || retinas(i, 3) < 2 %choose not to look at animals that have fewer than this number of each cell type
        continue
    end
    
    name = retinas(i, 1);
    fullname = ['SHOKR' num2str(name) 'L'];
    c_s{end + 1, 1} = fullname; 
    c_s{end, 2} = cell(0, 1);
    c_s{end, 3} = cell(0, 1);
    for j = 1:size(ups, 1)
        name_j = ups{j, 1};
        if contains(name_j, fullname)
            put = c_s{end, 2};
            put{end+1, 1} = ups{j, 2};
            c_s{end, 2} = put;
        end
    end
    
    for j = 1:size(downs, 1)
        name_j = downs{j, 1};
        if contains(name_j, fullname)
            put = c_s{end, 3};
            put{end+1, 1} = downs{j, 2};
            c_s{end, 3} = put;
        end
    end
end

% Now grab information about the firing rates of UP and DOWN cells on an animal by animal basis
windowTime = 25;
xvals = [windowTime/2:windowTime:3000-windowTime/2];
priorityOrder_1 = {'Flash'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; '_all'; 'Extracellular'}; %child field name should be one of these. Priority is in given order
c_iDiff = zeros(1, size(c_s, 1));
umf  = []; dmf = [];
for c = 1:size(c_s, 1)
   u = c_s{c, 2};
   d = c_s{c, 3};
   
   upPSTH = [];
   downPSTH = [];
   
   %up cells first
   for i = 1:size(u, 1)
       struct_i = u{i};
       
       analyses = fieldnames(struct_i.Analyses_Completed);
       analysisToUse_1 = findAnalysisToUse(priorityOrder_1, analyses);
       if ~strcmp(analysisToUse_1, 'NoneFound')
           analysisToUse_2 = findAnalysisToUse(priorityOrder_2, fieldnames(struct_i.Analyses_Completed.(analysisToUse_1)));
       else
           analysisToUse_2 = 'NoneFound';
       end
       
       
       if strcmp(analysisToUse_2, 'NoneFound')
           noAnalyses(end + 1) = i;
           continue
       end
       
       loc = struct_i.Analyses_Completed.(analysisToUse_1).(analysisToUse_2);
       
       preTime_i = loc.meta.preTime;
       stimTime_i = loc.meta.stimTime;
       tailTime_i = loc.meta.tailTime;
       
       totalTime = preTime_i + stimTime_i + tailTime_i;
       
       if totalTime ~= 3000
           continue
       end
       spikeTimes_i = loc.Analysis_Results.allSpikeTimes;
       
       numBins = numel(0:windowTime:(totalTime-windowTime));
       results_i = zeros(numel(spikeTimes_i), numBins);
       
       %average across all epochs
       for j = 1:numel(spikeTimes_i)
           s = spikeTimes_i{j};
           psth_j = zeros(1, numBins);
           for k = 0:numBins-1
               spikesInWindow = numel(s(s>k*windowTime & s<(k*windowTime+windowTime)))/(windowTime/1000);
               psth_j(k+1) = spikesInWindow;
           end
           results_i(j, :) = psth_j;
       end
       
       
       if j == 1
           meanResults = results_i;
       else
           meanResults = mean(results_i);
       end
       upPSTH = [upPSTH; meanResults];
   end
   
   %down cells second
   for i = 1:size(d, 1)
       struct_i = d{i};
       
       analyses = fieldnames(struct_i.Analyses_Completed);
       analysisToUse_1 = findAnalysisToUse(priorityOrder_1, analyses);
       if ~strcmp(analysisToUse_1, 'NoneFound')
           analysisToUse_2 = findAnalysisToUse(priorityOrder_2, fieldnames(struct_i.Analyses_Completed.(analysisToUse_1)));
       else
           analysisToUse_2 = 'NoneFound';
       end
       
       
       if strcmp(analysisToUse_2, 'NoneFound')
           noAnalyses(end + 1) = i;
           continue
       end
       
       loc = struct_i.Analyses_Completed.(analysisToUse_1).(analysisToUse_2);
       
       preTime_i = loc.meta.preTime;
       stimTime_i = loc.meta.stimTime;
       tailTime_i = loc.meta.tailTime;
       
       totalTime = preTime_i + stimTime_i + tailTime_i;
       
       if totalTime ~= 3000
           continue
       end
       spikeTimes_i = loc.Analysis_Results.allSpikeTimes;
       
       numBins = numel(0:windowTime:(totalTime-windowTime));
       results_i = zeros(numel(spikeTimes_i), numBins);
       
       %average across all epochs
       for j = 1:numel(spikeTimes_i)
           s = spikeTimes_i{j};
           psth_j = zeros(1, numBins);
           for k = 0:numBins-1
               spikesInWindow = numel(s(s>k*windowTime & s<(k*windowTime+windowTime)))/(windowTime/1000);
               psth_j(k+1) = spikesInWindow;
           end
           results_i(j, :) = psth_j;
       end
       
       
       if j == 1
           meanResults = results_i;
       else
           meanResults = mean(results_i);
       end
       downPSTH = [downPSTH; meanResults];
       
   end
   
   %Now, collect average data for this animal
   upFirst100ms = upPSTH(:, xvals>1000 & xvals < 1100);
   downFirst100ms = downPSTH(:, xvals>1000 & xvals<1100);
   
   maxUpFirst100ms = max(upFirst100ms, [], 2);
   maxDownFirst100ms = max(downFirst100ms, [],2);
   
   %c_iDiff stores the average difference between UP and DOWN cell firing
   %rates (first 100ms of flash response) on an animal by animal basis
   c_iDiff(c) = mean(maxDownFirst100ms) - mean(maxUpFirst100ms); %mean difference in first 100ms for cells of each type for each animal analyzed
   umf = [umf;maxUpFirst100ms]; %umf = max firing rate of up cells for all animals
   dmf = [dmf;maxDownFirst100ms]; %dmf = max firing rate of down cells for all animals
end

figure
title('Mean Firing Rate of UP Cells for Each Animal')
hold on
histogram(umf)
figure
title('Mean Firing Rate of DOWN Cells for Each Animal')
histogram(dmf)

figure
title('Difference in Max Firing Rate for UP and DOWN cells for Each Animal')
histogram(c_iDiff)