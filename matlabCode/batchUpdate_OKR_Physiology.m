load('C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\batch\cellList_OKR.mat')
analysisOfInterest = 'Moving_Bar_Speed_VarStimTime';
for i = 1:size(cellList)
    cell_i= cellList{i, 1};
    for j = 1:numel(cnames)
        cell_j = cnames{j};
        found = 0;
        switch cell_i
            case cell_j
                found = 1;
        end
        if found
            break
        end
    end
    if ~found
        continue
    end
    
    struct_i = cellList{i, 2};
    fpath = struct_i.fileLocation;
    fname = [fpath cell_i '_Analysis.mat'];
    load(fname) %loads CellParameters
    
    slashes = find(fname == '\');
    subfolder = fname(1:slashes(end-1));
    epochsFname = [subfolder, cell_i, '_ClarinetExport.mat'];
    load(epochsFname) %loads Epochs
    
    %find all instances of the analysis of interest
    analyses = CellParameters.Analyses_Completed.(analysisOfInterest);
    fs = fieldnames(analyses);
    for j = 1:numel(fs) %for each instance of the analysis of interest
        a_j = analyses.(fs{j}); %take the analysis substructure
        epochs_j = a_j.Analysis_Results.EpochNumbers; %Grab epoch numbers
        stimTimes = zeros(1, numel(epochs_j));
        for k = 1:numel(epochs_j)
            stimTimes(k) = epochs(epochs_j(k)).meta.currentStimTime; %grab relevant data from epochs
        end
        a_j.Analysis_Results.allStimTimes = stimTimes; %fill the analysis substructure
        analyses.(fs{j}) = a_j; %add the analysis sub structure back to the analyses 
    end
    
    CellParameters.Analyses_Completed.(analysisOfInterest) = analyses; %add analyses back to cell parameters
    save(fname, 'CellParameters')
    cellList{i, 2} = CellParameters; %syncs cellList with CellParameters
    cellList{i, 2}.updateHTML = 1;   
end
save('C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\batch\cellList_OKR.mat', 'cellList');
    