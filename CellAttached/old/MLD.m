function outputStruct = MLD(stimuliOrientations, cdfs, runsPerStimulus)
%Maximum Likelihood Decoder

%Inputs:
%stimuliOrientations - vector specifying stimulus values to be tested
%cdfs - 1xn cell array. Each cell contains a matrix of the cdfs of a cell to be used in the analysis. n = number of cells
%runsPerStimuli - integer number specifying how many times to test each
%stimulus orientation

%generate PDFs from each CDF (you'll use both pdfs and cdfs in the calculation)   
    
pdfs = cell(size(cdfs));
for c = 1:numel(cdfs)
    cell_c = cdfs{c};
    pdfs_c = zeros(size(cell_c));
    for o = 1:numel(stimuliOrientations)
        cdf_c = cell_c(o, :);
        pdf = diff([0 cdf_c]);
        pdf(end)= 1 - sum(pdf(1:end-1));
        pdfs_c(o, :) = pdf;
    end
    pdfs{c} = pdfs_c;
end


MLD_error = zeros(runsPerStimulus, 8); %store the errors by true stimulus angle
MLD_errorVertComponent = zeros(runsPerStimulus, 8); %store the error in the vertical component of the prediction by true stimulus angle
falseAlarms = cell(numel(stimuliOrientations), 1);
trueHits = zeros(numel(stimuliOrientations), 1);
falseAlarmCount = zeros(numel(stimuliOrientations), 1);
spikeResponses = zeros(numel(cdfs), numel(stimuliOrientations), runsPerStimulus);
trueStims = zeros(1, numel(stimuliOrientations)*runsPerStimulus);
predictedStims = zeros(1, numel(stimuliOrientations)*runsPerStimulus);

for i = 1:numel(stimuliOrientations)
    trueStim = stimuliOrientations(i); %grab the true stimulus
    trueStimVertComponent = sin(deg2rad(trueStim)); %compute the vertical component of the stimulus
        
    for c = 1:numel(cdfs)
        cdf = cdfs{c}; %grab the empirical CDFs for this stimulus and cell
        cdf = cdf(i, :);
        
        %generate n random responses for each cell, inserting correlated
        %noise if needed
        spikeResponses(c, i, :) = randNumSpikesGenerator(cdf, runsPerStimulus);   
    end
    
    for n = 1:runsPerStimulus
        %test the likelihood of every stimulus given the response
        likelihoods_n = zeros(1, numel(stimuliOrientations));
        
        for j = 1:numel(stimuliOrientations)
            %for each possible stimulus:
            probabilities = zeros(1, numel(cdfs)); %fill with probabilities of the given response for each cell involved
            
            %for each cell:
            for c = 1:numel(cdfs)
                p_c = pdfs{c}; %grab the PDFs for this cell
                pdf = p_c(j, :); %grab the pdf for this stimulus
                response_c = spikeResponses(c, i, n); %grab the response to this stimulus
               try
                pResponse = pdf(round(response_c*(numel(pdf)))); %find the probability of this response given the test stimulus
               catch
               end
                probabilities(c) = pResponse;
            end
            
            %after finding the probability of the test stimulus for each
            %cell's response, calculate the likelihood of the test stimulus
            likelihood = prod(probabilities);
            likelihoods_n(j) = likelihood;
        end
        
        %after calculating the likelihood of each stimulus, find the one
        %with the highest likelihood... this is your prediction
        [~, maxLikelihoodIndx] = max(likelihoods_n);
        predictedStim = stimuliOrientations(maxLikelihoodIndx);
        predictedVertComponent = sin(deg2rad(predictedStim));
        
        trueStims((i-1)*runsPerStimulus+n) = trueStim;
        predictedStims((i-1)*runsPerStimulus+n) = predictedStim;
        
        %compare with the true stimulus and count false alarms and correct
        %hits
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

%confusion matrix stores the frequency with which the decoder guessed
%each of the possible stimuli, given the true stimulus
confusionMatrix = confusionmat(trueStims,predictedStims);

accuracy = sum(trueStims == predictedStims)/numel(trueStims);

%store everything in an output matrix and return
outputStruct = struct();
outputStruct.accuracy = accuracy;
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
outputStruct.confusionMatrix = confusionMatrix;
outputStruct.stimuli = stimuliOrientations;
outputStruct.runsPerStimuli = runsPerStimulus;
outputStruct.simulatedResponses = spikeResponses;
outputStruct.trueStims = trueStims;
outputStruct.predictedStims = predictedStims;
end

