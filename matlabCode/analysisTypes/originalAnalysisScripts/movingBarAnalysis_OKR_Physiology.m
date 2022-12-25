function [processedDataStruct, analysisFigure] = movingBarAnalysis_OKR_Physiology(data,...
    recordingType, analysisType, epochsSelected, cellID)
%% Analysis script for moving bar stim
    processedDataStruct = struct();

    %% Extracellular (spikes) analysis
    if strcmp(recordingType, 'Extracellular')

        %create a Nx2 cell array to hold responses oreintation by orientation
        %The first column will hold a single number that corresponds to the orientation
        %The second column will hold a 1xM array that contains the
        %responses at that oreintation. N corresponds to the number of
        %orientations (usually 8), and M corresponds to the number of
        %epochs at each orientation (often 5).
        allOrientations = data.(cellID).epochs(epochsSelected(1)).meta.orientation; %in degrees, 1xN vector of all orientations used
        allOrientations = convertToRetinaAngle(allOrientations, data.(cellID).fileLocation); %convert angles based on how projector is flipped
        responsesByOrientation = cell(numel(allOrientations), 2); %initialize the Nx2 cell array        
        responsesByOrientation(:, 1) = num2cell(allOrientations');%fill the first column of the cell array with all orientations (in degrees)
        
        speed = data.(cellID).epochs(epochsSelected(1)).meta.speed; %degrees/second        
        barWidth = data.(cellID).epochs(epochsSelected(1)).meta.barWidth;
        barHeight = data.(cellID).epochs(epochsSelected(1)).meta.barHeight;
        
        processedDataStruct.Moving_Bar.Extracellular.allSpikeTimes = {}; %contains a row vector for each epoch which shows the time (in ms) of each spike occurance during that epoch    
        
        orientationByEpoch = []; %record which epoch had which orientation.
        %Loop through the selected epochs one at a time.
        for i = 1:numel(epochsSelected)
            branch_i = data.(cellID).epochs(epochsSelected(i));
            
            %Pull out the recording trace
            trace = branch_i.epoch;
            
            %get relevant parameters
            currentOrientation = branch_i.meta.currentOrientation;%deg
            currentOrientation = convertToRetinaAngle(currentOrientation, data.(cellID).fileLocation); %convert angles based on how projector is flipped
            preTime = branch_i.meta.preTime; %ms
            stimTime = branch_i.meta.stimTime; %ms
            sampleRate = branch_i.meta.sampleRate; %observations per second
               
            orientationByEpoch = [orientationByEpoch, currentOrientation];
            
            %get spike times
            spikeData = SpikeDetector_SK(trace);
            %spikeData = simpleSpikeDetector_SH(trace, sampleRate); %use when original spike detector doesn't work

            spikeTimes_SR = spikeData.sp; %times in units of sampleRate;
            spikeTimes = (spikeTimes_SR./sampleRate)*1000; %convert to ms;
            
            processedDataStruct.Moving_Bar.Extracellular.allSpikeTimes{i} = spikeTimes;

            %count spikes that happen during the stim
            stimSpikeIndx = find(spikeTimes > preTime &...
                spikeTimes < (preTime + stimTime));
            numSpikes = numel(stimSpikeIndx);
            
            %append numSpikes to relevant position in
            %responsesByOrientation cell array
            indexToUse = find(cell2mat(responsesByOrientation(:,1)) == currentOrientation);
            
            if isempty(indexToUse) %if the current orienatation doesn't exist, add it to the bottom
                responsesByOrientation{end + 1, 1} = currentOrientation;
                responsesByOrientation{end, 2} = [];
                indexToUse = find(cell2mat(responsesByOrientation(:,1)) == currentOrientation);
                allOrientations = [allOrientations, currentOrientation];
            end
            
            responsesByOrientation{indexToUse, 2} = [responsesByOrientation{indexToUse, 2}, numSpikes];
        end
        
        %fill these fields once... they'll take the value from the last
        %epoch but really it should be the same across epochs
        processedDataStruct.Moving_Bar.Extracellular.preTime = preTime; %time before stimulus onset
        processedDataStruct.Moving_Bar.Extracellular.stimTime = stimTime; %time of stimulus
        processedDataStruct.Moving_Bar.Extracellular.tailTime = branch_i.meta.tailTime; %time after stimulus offset
        processedDataStruct.Moving_Bar.Extracellular.sampleRate = sampleRate; %observations per second
        processedDataStruct.Moving_Bar.Extracellular.barWidth = barWidth;
        processedDataStruct.Moving_Bar.Extracellular.barHeight = barHeight;
        processedDataStruct.Moving_Bar.Extracellular.barIntensity = branch_i.meta.barIntensity; %amplitude of the grating (normalized by maximum possible based on meanIntensity)
        processedDataStruct.Moving_Bar.Extracellular.backgroundIntensity = branch_i.meta.backgroundIntensity; %background intensity that occurs during pre and tail time
        processedDataStruct.Moving_Bar.Extracellular.centerOffset = branch_i.meta.centerOffset;
        processedDataStruct.Moving_Bar.Extracellular.orientationByEpoch = orientationByEpoch;

        %calculate mean and standard deviation spike number for each direction
        meanByOrientation = [];
        stdByOrientation = [];
        for i = 1:numel(allOrientations)
            meanByOrientation = [meanByOrientation mean(responsesByOrientation{i,2})];
            stdByOrientation = [stdByOrientation std(responsesByOrientation{i, 2})];
        end
        
        %calculate mean vector
        allOrientationsRad = deg2rad(allOrientations);
        [x,y] = pol2cart(allOrientationsRad, meanByOrientation);
        xSUM = sum(x);
        ySUM = sum(y);
        [meanTheta, meanRho] = cart2pol(xSUM, ySUM);
        meanTheta = rad2deg(meanTheta);
        
        %calculate DSI... DSI = length of preferred vector/sum(response in
        %all directions).
        totalSum = sum(meanByOrientation); %sum from all directions
        DSI = meanRho/totalSum; %vector length over sum of all the directions
        
        %add to processed data struct
        processedDataStruct.Moving_Bar.Extracellular.EpochNumbers = epochsSelected;
        processedDataStruct.Moving_Bar.Extracellular.Orientation = allOrientations;
        processedDataStruct.Moving_Bar.Extracellular.Speed = speed;
        processedDataStruct.Moving_Bar.Extracellular.spikesByOrientation = responsesByOrientation;
        processedDataStruct.Moving_Bar.Extracellular.mean_spikesByOrientation = meanByOrientation;
        processedDataStruct.Moving_Bar.Extracellular.std_spikesByOrientation = stdByOrientation;
        processedDataStruct.Moving_Bar.Extracellular.PreferredDirection = meanTheta;
        processedDataStruct.Moving_Bar.Extracellular.MeanVectorLength = meanRho;
        processedDataStruct.Moving_Bar.Extracellular.DSI = DSI;
        
        %add metadata from last epoch to the structure, but remove
        %irrelivent fields.
        irreliventFields = {'bathTemperature', 'epochNum','currentOrientation'};
        metaToAdd = branch_i.meta;
        editedMeta = rmfield(metaToAdd, irreliventFields);
        processedDataStruct.Moving_Bar.Extracellular.meta = editedMeta;
        
        %% create analysis figure
        analysisFigure = figure;
        polarwitherrorbar(deg2rad(allOrientations), meanByOrientation, stdByOrientation);
        hold on
        polarplot([0 deg2rad(meanTheta)],[0 meanRho], '-k', 'LineWidth', 2);
        title(strcat(cellID, ' Moving Bar'));
        caption = ['Speed = ' num2str(round(speed, 2))...
            '; Preferred Angle = ' num2str(round(meanTheta, 1)) '; Vector Length = ' num2str(round(meanRho, 2))...
            '; DSI = ' num2str(round(DSI, 3))];
        a = annotation('textbox',[0.02 0.06 0 0], 'String',caption,'FitBoxToText','on');
        a.FontSize = 9;
        hold off
        
    end
end