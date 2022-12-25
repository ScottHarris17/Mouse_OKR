%use this to plot individual epochs from single unit recordings in the OKR
%project without having to use data viewer. Give cell name and epoch number
%as input values. Scroll to "plot" section to edit figure.

cellID = 'SH15Lc2'; %choose the cell
epochNums = 22; %choose the epoch numbers to plot

cellListFName = "C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/batch/cellList_OKR.mat";
%% load data
load(cellListFName);

for i = 1:size(cellList, 1)
    if strcmp(cellList{i, 1}, cellID)
        cellDataStruct = cellList{i,2};
        break
    end
end

for i = 1:numel(epochNums)
    epochNumber = epochNums(i);
    epochs = getClarinetExport(cellDataStruct);

    meta = epochs(epochNumber).meta;
    trace = epochs(epochNumber).epoch;
    trace = trace - mean(trace);
    totalTime = meta.preTime + meta.stimTime + meta.tailTime;
    time = linspace(0, totalTime, numel(trace));


%% plot
    figure
    hold on
    plot(time, trace, 'Color', [190 200 120]/255)
    set(gca,'visible','off')
    title(['Tonic Grating Response'])
    plot([meta.preTime, meta.preTime+meta.stimTime],...
        [mean(trace) + 12*std(trace), mean(trace) + 12*std(trace)],...
        'Color', [0 0 0]/255, 'LineWidth', 2)
    xlabel('ms')
    ylabel('mV')
    hold off
end

%% get plots for a specific current orientation
% angle = 270;
% for m = 1:size(epochs, 2)
%     try
%         if epochs(m).meta.currentOrientation == angle...
%                 && strfind(epochs(m).meta.displayName
%             meta = epochs(m).meta;
%             trace = epochs(m).epoch;
%             trace = trace - mean(trace);
%             totalTime = meta.preTime + meta.stimTime + meta.tailTime;
%             time = linspace(0, totalTime, numel(trace));
%             figure
%             plot(time, trace, 'k')
%             hold on
%             plot([meta.preTime, meta.preTime+meta.stimTime],...
%             [mean(trace) + 12*std(trace), mean(trace) + 12*std(trace)],...
%             'Color', [222 53 118]/255, 'LineWidth', 2)
%             set(gca,'visible','off')
%             hold off
%         end
%     catch
%         continue
%     end
% end
% 
