function [processedDataStruct, analysisFigure] = movingGratingSpeedAnalysis_OKR_Physiology(data,...
    recordingType, analysisType, epochsSelected, cellID)
    %% Analysis script for moving grating speed stim
    processedDataStruct = struct();
    
     %% Extracellular (spikes) analysis
    if strcmp(recordingType, 'Extracellular')
        
        %create a Nx2 cell array to hold responses speed by speed
        %The first column will hold a single number that corresponds to the
        %speed. The second column will hold a 1xM array that contains the
        %responses at each speed. N corresponds to the number of
        %speeds, and M corresponds to the number of
        %epochs at each speed (usually 5).
        allSpeeds = data.(cellID).epochs(epochsSelected(1)).meta.speed; %in degrees, 1xN vector of all speeds used
        responsesBySpeed = cell(numel(allSpeeds), 2); %initialize the Nx2 cell array        
        responsesBySpeed(:, 1) = num2cell(allSpeeds');%fill the first column of the cell array with all speeds (degrees/second)
        
        orientation = data.(cellID).epochs(epochsSelected(1)).meta.orientation; %degrees
        orientation = convertToRetinaAngle(orientation, data.(cellID).fileLocation); %convert angles based on how projector is flipped
        spatialFrequency = data.(cellID).epochs(epochsSelected(1)).meta.spatialFrequency;%degrees/cycle
        
        processedDataStruct.Moving_Grating_Speed.Extracellular.allSpikeTimes = {}; %contains a row vector for each epoch which shows the time (in ms) of each spike occurance during that epoch    
        
        speedByEpoch = [];
        
        %Loop through the selected epochs one at a time.
        for i = 1:numel(epochsSelected)
            branch_i = data.(cellID).epochs(epochsSelected(i));
            
            %Pull out the recording trace
            trace = branch_i.epoch;
            
            %get relevant parameters
            currentSpeed = branch_i.meta.currentSpeed;%deg/s
            preTime = branch_i.meta.preTime; %ms
            stimTime = branch_i.meta.stimTime; %ms
            sampleRate = branch_i.meta.sampleRate; %observations per second
               
            speedByEpoch = [speedByEpoch, currentSpeed];
            
            %get spike times
            spikeData = SpikeDetector_SK(trace);
            spikeTimes_SR = spikeData.sp; %times in units of sampleRate;
            spikeTimes = (spikeTimes_SR./sampleRate)*1000; %convert to ms;
            
            processedDataStruct.Moving_Grating_SpatialFrequency.Extracellular.allSpikeTimes{i} = spikeTimes;

            %count spikes that happen during the stim
            stimSpikeIndx = find(spikeTimes > preTime &...
                spikeTimes < (preTime + stimTime));
            numSpikes = numel(stimSpikeIndx);
            
            %append numSpikes to relevant position in
            %responsesByOrientation cell array
            indexToUse = find(cell2mat(responsesBySpeed(:,1)) == currentSpeed);
           
            if isempty(indexToUse) %if the current speed doesn't exist, add it to the bottom
                responsesBySpeed{end + 1, 1} = currentSpeed;
                responsesBySpeed{end, 2} = [];
                indexToUse = find(cell2mat(responsesBySpeed(:,1)) == currentSpeed);
                allSpeeds = [allSpeeds, currentSpeed];
            end
            
            responsesBySpeed{indexToUse, 2} = [responsesBySpeed{indexToUse, 2}, numSpikes];
        end

        processedDataStruct.Moving_Grating_Speed.Extracellular.preTime = preTime; %time before stimulus onset
        processedDataStruct.Moving_Grating_Speed.Extracellular.stimTime = stimTime; %time of stimulus
        processedDataStruct.Moving_Grating_Speed.Extracellular.tailTime = branch_i.meta.tailTime; %time after stimulus offset
        processedDataStruct.Moving_Grating_Speed.Extracellular.sampleRate = sampleRate; %observations per second
        processedDataStruct.Moving_Grating_Speed.Extracellular.meanIntensity = branch_i.meta.meanIntensity; %mean intensity of the grating
        processedDataStruct.Moving_Grating_Speed.Extracellular.amplitude = branch_i.meta.amplitude; %amplitude of the grating (normalized by maximum possible based on meanIntensity)
        processedDataStruct.Moving_Grating_Speed.Extracellular.backgroundIntensity = branch_i.meta.backgroundIntensity; %background intensity that occurs during pre and tail time
        processedDataStruct.Moving_Grating_Speed.Extracellular.speedByEpoch = speedByEpoch;

        %calculate mean and standard deviation spike number for each speed
        meanBySpeed = [];
        stdBySpeed = [];
        for i = 1:numel(allSpeeds)
            meanBySpeed = [meanBySpeed mean(responsesBySpeed{i,2})];
            stdBySpeed = [stdBySpeed std(responsesBySpeed{i, 2})];
        end
        

        
        %add to processed data struct
        processedDataStruct.Moving_Grating_Speed.Extracellular.EpochNumbers = epochsSelected;
        processedDataStruct.Moving_Grating_Speed.Extracellular.Speed = allSpeeds;
        processedDataStruct.Moving_Grating_Speed.Extracellular.Orientation = orientation;
        processedDataStruct.Moving_Grating_Speed.Extracellular.SpatialFrequency = spatialFrequency;
        processedDataStruct.Moving_Grating_Speed.Extracellular.spikesByspeed = responsesBySpeed;
        processedDataStruct.Moving_Grating_Speed.Extracellular.mean_spikesBySpeed = meanBySpeed;
        processedDataStruct.Moving_Grating_Speed.Extracellular.std_spikesBySpeed = stdBySpeed;
        
        %add metadata from last epoch to the structure, but remove
        %irrelivent fields.
        irreliventFields = {'bathTemperature', 'epochNum','currentSpeed'};
        metaToAdd = branch_i.meta;
        editedMeta = rmfield(metaToAdd, irreliventFields);
        processedDataStruct.Moving_Grating_Speed.Extracellular.meta = editedMeta;
    end
    %% Build analysis figure
    X = allSpeeds;
    Y = meanBySpeed;
    errorTop = Y + stdBySpeed;
    errorBottom = Y - stdBySpeed;

    analysisFigure = figure();
    bar(X, Y)
    hold on
    erPlot = errorbar(X, Y, errorBottom, errorTop);
    erPlot.LineStyle = 'none';
    title(strcat(cellID, 'Moving Grating Speed'));
    xlabel('speed (deg/s)')
    ylabel('Mean Response (# spikes)');
    caption = ['Orientation = ' num2str(orientation) '; Spatial Frequency = ' num2str(round(spatialFrequency, 2))];
    a = annotation('textbox',[0.5 0.65 0 0], 'String',caption,'FitBoxToText','on');
    a.FontSize = 9;
    hold off
end
