function   [processedDataStruct, analysisFigure] = FlashAnalysis_OKR_Physiology(data, recordingType, epochsSelected, cellID, cellType)
    %This function processes data from a Flash stimulus. It is called by
    %analyzer_OKR_Physiology.m if the user has chosen to analyze Flash
    %epochs from the dataViewer gui. Unlike other analysis types, Flash 
    %analysis can be run on epochs that used a variety of stimuli if the 
    %user desires, though this is not necessarily advised.
    %INPUTS:
    %data = the full original data structure with all cells' data in it
    %recordingType = a string representing the recording type (e.g.
    %'Extracellular')
    %epochsSelected = string containing all the epoch numbers that the user
    %would like to anlayze (e.g. '1 2 3 10 34')
    %cellID = a string representing the name of the cell that the user is
    %analyzing
    %OUTPUTS:
    %processedDataStruct = a structure containing the processed data
    
    %% Modifiable Parameters
    transientThreshold = 100; %in ms, threshold cutoff between transient and sustained flash responses

    processedDataStruct = struct();
    %% Extracellular (spikes) analysis
    if strcmp(recordingType, 'Extracellular')
        
        %initialize processed data structure;
        processedDataStruct.Flash.Extracellular.EpochNumbers = epochsSelected;%epochs analyzed in this analysis run
        processedDataStruct.Flash.Extracellular.transientThreshold = transientThreshold; %copy of the transient threshold set by user above
        processedDataStruct.Flash.Extracellular.allSpikeTimes = {}; %contains a row vector for each epoch which shows the time (in ms) of each spike occurance during that epoch
        processedDataStruct.Flash.Extracellular.baseline_firingRates = []; %list of firing rates from all epochs during baseline
        processedDataStruct.Flash.Extracellular.ON_firingRates = []; %list of firing rates from all epochs during stim on (during the flash)
        processedDataStruct.Flash.Extracellular.ON_transient_firingRates = []; %list of firing rates from all epochs between stim on and transient threshold (set above)
        processedDataStruct.Flash.Extracellular.ON_sustained_firingRates = []; %list of firing rates from all epochs during stim on and after transient threshold
        processedDataStruct.Flash.Extracellular.OFF_firingRates = []; %list of firing rates from all epochs during tail time
        processedDataStruct.Flash.Extracellular.allSpikeTimes = {};
        
        %Loop through the selected epochs one at a time.
        for i = 1:numel(epochsSelected)
            branch_i = data.(cellID).epochs(epochsSelected(i));
            
            %Pull out the recording trace
            trace = branch_i.epoch;

            %Pull out the relevant metadata for analysis
            preTime = branch_i.meta.preTime; %in ms
            stimTime = branch_i.meta.stimTime; %in ms
            tailTime = branch_i.meta.tailTime; %in ms
            sampleRate = branch_i.meta.sampleRate; %observations per second
            
            %Extract spikes using SK's algorithm
            spikeData = SpikeDetector_SK(trace);
            %spikeData = simpleSpikeDetector_SH(trace, sampleRate); %when original doesn't work
            
            spikeTimes_SR = spikeData.sp; %times in units of sampleRate;
            
            spikeTimes = (spikeTimes_SR./sampleRate)*1000; %convert to ms;

            
            processedDataStruct.Flash.Extracellular.allSpikeTimes{i} = spikeTimes; %add to struct
            
            %find spike times that happened during the pretime
            preSpikeIndx = find(spikeTimes < preTime);
            stimSpikeIndx = find(spikeTimes > preTime &...
                spikeTimes < (preTime + stimTime));
            tailSpikeIndx = find(spikeTimes > (preTime + stimTime));
            
            %calculate firing rates over each part of the epoch
            baselineFR = numel(preSpikeIndx)/(preTime/1000); %firing rate in hz;
            ON_FR = numel(stimSpikeIndx)/(stimTime/1000);
            OFF_FR = numel(tailSpikeIndx)/(tailTime/1000);
            
            %Get a bit more detailed by looking for ON transient/sustained
            ONtransientSpikes = spikeTimes(stimSpikeIndx) < (stimTime + transientThreshold);
            ONsustainedSpikes = spikeTimes(stimSpikeIndx) > (stimTime + transientThreshold);
            
            numONtransient = numel(find(ONtransientSpikes));
            numONsustained = numel(find(ONsustainedSpikes));
            
            ONtransientFR = numONtransient/(transientThreshold/1000);
            ONsustainedFR = numONsustained/((stimTime - transientThreshold)/1000);
            
            %add to processedDataStruct
            processedDataStruct.Flash.Extracellular.baseline_firingRates = [processedDataStruct.Flash.Extracellular.baseline_firingRates; baselineFR];
            processedDataStruct.Flash.Extracellular.ON_firingRates = [processedDataStruct.Flash.Extracellular.ON_firingRates; ON_FR];
            processedDataStruct.Flash.Extracellular.ON_transient_firingRates = [processedDataStruct.Flash.Extracellular.ON_transient_firingRates; ONtransientFR];
            processedDataStruct.Flash.Extracellular.ON_sustained_firingRates = [processedDataStruct.Flash.Extracellular.ON_sustained_firingRates; ONsustainedFR];
            processedDataStruct.Flash.Extracellular.OFF_firingRates = [processedDataStruct.Flash.Extracellular.OFF_firingRates; OFF_FR];
        end
        
      
        %% calculate means and standard deviations for each metric
        avg_baseline = mean(processedDataStruct.Flash.Extracellular.baseline_firingRates);
        std_baseline = std(processedDataStruct.Flash.Extracellular.baseline_firingRates);
        avg_ON_FR = mean(processedDataStruct.Flash.Extracellular.ON_firingRates);
        std_ON_FR = std(processedDataStruct.Flash.Extracellular.ON_firingRates);
        avg_ONtransient_FR = mean(processedDataStruct.Flash.Extracellular.ON_transient_firingRates);
        std_ONtransient_FR = std(processedDataStruct.Flash.Extracellular.ON_transient_firingRates);
        avg_ONsustained_FR = mean(processedDataStruct.Flash.Extracellular.ON_sustained_firingRates);
        std_ONsustained_FR = std(processedDataStruct.Flash.Extracellular.ON_sustained_firingRates);
        avg_OFF_FR = mean(processedDataStruct.Flash.Extracellular.OFF_firingRates);
        std_OFF_FR = std(processedDataStruct.Flash.Extracellular.ON_firingRates);
        
        %and then add all of these back into the processedDataStruct     
        processedDataStruct.Flash.Extracellular.average_baseline_firingRate = avg_baseline;
        processedDataStruct.Flash.Extracellular.std_baseline_firingRate = std_baseline;
        processedDataStruct.Flash.Extracellular.average_ON_firingRate = avg_ON_FR;
        processedDataStruct.Flash.Extracellular.std_ON_firingRate = std_ON_FR;
        processedDataStruct.Flash.Extracellular.average_ON_transient_firingRate = avg_ONtransient_FR;
        processedDataStruct.Flash.Extracellular.std_ON_transient_firingRate = std_ONtransient_FR;
        processedDataStruct.Flash.Extracellular.average_ON_sustained_firingRate = avg_ONsustained_FR;
        processedDataStruct.Flash.Extracellular.std_ON_sustained_firingRate = std_ONsustained_FR;
        processedDataStruct.Flash.Extracellular.average_OFF_firingRate = avg_OFF_FR;
        processedDataStruct.Flash.Extracellular.std_OFF_firingRate = std_OFF_FR;
        
        %add metadata from last epoch to the structure, but remove
        %irrelivent fields.
        irreliventFields = {'bathTemperature'};
        metaToAdd = branch_i.meta;
        editedMeta = rmfield(metaToAdd, irreliventFields);
        processedDataStruct.Flash.Extracellular.meta = editedMeta;
        
        
        %% Create analysis figure. It is a bar plot with the average spike
        %rates from prestim, stim (on), and post stim (off) time points
        analysisFigure = figure;
        
        X = categorical({'baseline','ON','OFF'});
        X = reordercats(X,{'baseline','ON','OFF'});       
        Y = [avg_baseline, avg_ON_FR, avg_OFF_FR];
        bar(X, Y);
        hold on
       
        err = [std_baseline,std_ON_FR, std_OFF_FR];
        errorTop = err + Y;
        errorBottom = Y - err;
        
        erPlot = errorbar(X, Y, errorBottom, errorTop);
        erPlot.LineStyle = 'none';
        
        ylabel('Average Firing Rate');
        title(strcat(cellID, ' Flash - ', cellType));
        hold off
        
    end
end
    
