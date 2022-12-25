function outputStruct = maximumLikelihoodDecoder_old(stimuliOrientations, cdfsMat1, cdfsMat2, runsPerStimulus)
%Maximum Likelihood Decoder
MLD_error = zeros(runsPerStimulus, 8); %store the errors by true stimulus angle
MLD_errorVertComponent = zeros(runsPerStimulus, 8); %store the error in the vertical component of the prediction by true stimulus angle
falseAlarms = cell(8, 1);
trueHits = zeros(8, 1);
falseAlarmCount = zeros(8, 1);

for i = 1:numel(stimuliOrientations)
    trueStim = stimuliOrientations(i); %grab the true stimulus
    trueStimVertComponent = sin(deg2rad(trueStim)); %compute the vertical component of the stimulus
    cdf_i1 = cdfsMat1(i, :); %grab the empirical CDFs for this stimulus
    cdf_i2 = cdfsMat2(i, :);
    
    seeds1 = rand(1, runsPerStimulus); %compute random responses for all runs
    seeds2 = rand(1, runsPerStimulus);
    
    for n = 1:runsPerStimulus %adjust the random responses to the xFill value that is closes to them
        [~, indx1] = min((cdf_i1 - seeds1(n)).^2);
        [~, indx2] = min((cdf_i2 - seeds2(n)).^2);
                
        %test the likelihood of every stimulus given the response
        likelihoods = zeros(1, numel(stimuliOrientations));
        for j = 1:numel(stimuliOrientations)
            %for each possible stimulus:
            
            %calculate the PDF for this stimulus for the first group of
            %cells
            askStimCDF_1 = cdfsMat1(j, :);
            askStimPDF_1 = diff([0 askStimCDF_1]);
            askStimPDF_1(end) = 1 - sum(askStimPDF_1(1:end-1));
            
            %calculate the probability that this test stimulus elicited the
            %group 1 response
            p_askStim1 = askStimPDF_1(indx1);
            
            %grab the CDF for group 2 and compute the PDF
            askStimCDF_2 = cdfsMat2(j, :);
            askStimPDF_2 = diff([0 askStimCDF_2]);
            askStimPDF_2(end) = 1 - sum(askStimPDF_2(1:end-1));
            
            %calculate the probability that this test stimulus elicited the
            %group 2 response
            p_askStim2 = askStimPDF_2(indx2);
            
            %calculate the likelihood, which is the product of the two
            %probabilities (it assumes independence, i.e. no noise
            %correlation).
            likelihood = p_askStim1 * p_askStim2;
            likelihoods(j) = likelihood;
        end
        
        %for each mock response, find the angle with stimulus with the
        %highest likelihood of eliciting it
        [~, maxLikelihoodIndx] = max(likelihoods);
        predictedStim = stimuliOrientations(maxLikelihoodIndx);
        predictedVertComponent = sin(deg2rad(predictedStim));
        
        if predictedStim ~= trueStim
            falseAlarmCount(stimuliOrientations == predictedStim) = falseAlarmCount(stimuliOrientations == predictedStim) + 1;
            falseAlarms{stimuliOrientations == predictedStim} =  [falseAlarms{stimuliOrientations == predictedStim}, trueStim];
        else
            trueHits(i) = trueHits(i) + 1;
        end  
        
        %calculate the errors and store in the larger matrices
        error = abs(trueStim - predictedStim);
        if error > 180 %always get the smaller angle
            error = min([abs(360-trueStim - predictedStim), abs(trueStim - (360-predictedStim))]);
        end
        errorVertComponent = abs(trueStimVertComponent - predictedVertComponent);
        MLD_error(n, i) = error;
        MLD_errorVertComponent(n, i) = errorVertComponent;
    end
end

classificationMatrix = zeros(numel(stimuliOrientations));
for i = 1:size(falseAlarms)
    alarms_i = falseAlarms{i};
    for j = 1:numel(stimuliOrientations)
        orientation_j = stimuliOrientations(j);
        classificationMatrix(i, j) = sum(alarms_i == orientation_j)./runsPerStimulus;
        if i == j
            classificationMatrix(i, j) = trueHits(i)/runsPerStimulus;
        end
    end
end

outputStruct = struct();
outputStruct.errorByStimulus = MLD_error;
outputStruct.meanError = mean(MLD_error);
outputStruct.verticalComponentErrorByStimulus = MLD_errorVertComponent;
outputStruct.meanVerticalComponentError = mean(MLD_errorVertComponent);
outputStruct.trueHitCountByStimulus = trueHits;
outputStruct.trueHitRate = trueHits./runsPerStimulus;
outputStruct.falseAlarms = falseAlarms;
outputStruct.falseAlarmCountByStimulus = falseAlarmCount;
outputStruct.falseAlarmRate = falseAlarmCount./runsPerStimulus;
outputStruct.hitToMissRatio = trueHits./falseAlarmCount;
outputStruct.classificationMatrix = classificationMatrix;
outputStruct.stimuli = stimuliOrientations;
outputStruct.runsPerStimuli = runsPerStimulus;
end

