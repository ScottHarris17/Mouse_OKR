%Run this script first when doing current clamp tuning curves for many
%contrasts. After this use the extractVms.mlapp to extract underlying Vms
%of the traces

%This is the preprocessing script. It simply extracts the clarinet export
%data that corresponds to the tuning curves at different contrasts and
%saves it to a new folder in the same directory as this script.

ups = grabFromFilter('UP Cells');
downs = grabFromFilter('DOWN Cells');

allCells = [ups;downs];
%choose which protocols to look at
priorityOrder_1 = {'Moving_Bar'}; %first field name shosuld be one of these. Priority is in given order
priorityOrder_2 = {'_SP10'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'CurrentClamp', 'PotassiumSpikes', '_ff'};
mustHave_2b = {'CurrentClamp', 'PotassiumSpikes', '_Intensity0d1'};
cantHave_1 = [];
cantHave_2 = {'Cesium'}; %avoid anything with these tags
rigMandates = []; %no rig mandates

dataLocation = 'C:\Users\mrsco\Box\Completed_Analysis\WholeCell\ManyContrasts';
upSpikeCurves = [];
upVmCurves = [];
upSpikeDSIs = [];
upVmDSIs = [];
upSpikePreferredDirections = [];
upVmPreferredDirections = [];
upNames = {};
for i = 1:size(allCells)
    struct_i = allCells{i, 2};
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    [loc_lowC, pass_LowC] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2b, cantHave_2, rigMandates);
    
    if pass || pass_LowC
        continue
    end
    
    cname = struct_i.cellID;
    folderName = [dataLocation, '\', cname];
    
    if ~exist(folderName, 'dir')
        cellDataPath = struct_i.fileLocation;
        slashes = strfind(cellDataPath, '\');
        clarinetExportPath = [cellDataPath(1:slashes(end-1)) cname '_ClarinetExport.mat'];
        load(clarinetExportPath) %loads struct named epochs
        
        ffEpochNums = loc.Analysis_Results.EpochNumbers;
        Intensity0d1EpochNums = loc_lowC.Analysis_Results.EpochNumbers;
        
        ffEpochs = zeros(numel(ffEpochNums), numel(epochs(ffEpochNums(1)).epoch));
        for j = 1:numel(ffEpochNums)
            ffEpochs(j, :) = epochs(ffEpochNums(j)).epoch;
        end
        
        Intensity0d1Epochs = zeros(numel(Intensity0d1EpochNums), numel(epochs(Intensity0d1EpochNums(1)).epoch));
        for j = 1:numel(Intensity0d1EpochNums)
            Intensity0d1Epochs(j, :) = epochs(Intensity0d1EpochNums(j)).epoch;
        end
        
        ffEpochs = ffEpochs';
        Intensity0d1Epochs = Intensity0d1Epochs';
        
        mkdir(folderName)
        save([folderName, '\ffEpochs.mat'], 'ffEpochs')
        save([folderName, '\Intensity0d1Epochs.mat'], 'Intensity0d1Epochs')      
    end
end