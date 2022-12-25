function   [processedDataStruct, analysisFigure] = analyzer_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID, cellType)
    %link from dataViewerViewerApp_OKR to all analysis scripts
    %
    %INPUTS:
    %Data = the full original data structure with all cells' data in it
    %recordingType = a string representing the recording type (e.g.
    %'extracellular')
    %analysisType = a string representing the type of analysis to do (e.g.
    %'Moving Grating - Direction')
    %epochsSelected = a row vector containing all the epoch numbers that the user
    %would like to anlayze (e.g. [1 2 3 10 34])
    %cellID = a string representing the name of the cell that the user is
    %analyzing
    
    %OUTPUTS:
    %processedDataStruct = a structure containing the processed data
    
    disp(['Running ' analysisType ', ' recordingType ' from ' cellID])
    %% Call the desired analysisType function
    
        if strcmp(analysisType, 'Flash')
            [processedDataStruct, analysisFigure] = flashAnalysis_OKR_Physiology(data, recordingType, epochsSelected, cellID, cellType);
           
        elseif strcmp(analysisType, 'Bar Noise')
            [processedDataStruct, analysisFigure] = barNoiseAnalysis_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);
        
        elseif strcmp(analysisType, 'Moving Bar')
            [processedDataStruct, analysisFigure] = movingBarAnalysis_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);
            
        elseif strcmp(analysisType, 'Moving Grating Direction')  
            [processedDataStruct, analysisFigure] = movingGratingDirectionAnalysis_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);
        
        elseif strcmp(analysisType, 'Moving Grating Speed')
            [processedDataStruct, analysisFigure] = movingGratingSpeedAnalysis_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);

        elseif strcmp(analysisType, 'Moving Grating Spatial Frequency')  
            [processedDataStruct, analysisFigure] = movingGratingSpatialFrequencyAnalysis_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);           
       
        elseif strcmp(analysisType, 'Sinusosoidal Oscillation')
            [processedDataStruct, analysisFigure] = sinusoidalOscillation_OKR_Physiology(data, recordingType, analysisType, epochsSelected, cellID);
        
        end
    
        %% Add anything else you need to processedDataStruct
        processedDataStruct.cellType = cellType;
end
