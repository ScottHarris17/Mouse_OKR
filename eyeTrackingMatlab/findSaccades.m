function saccadeDetection = findSaccades(eyeTrace, frameRate)
%This function calculates the location of saccades in an eye trace
%inputs:
%     - eyeTrace - 1xn or nx1 vector of eye position across time
%     - frameRate - frame rate of the eye trace (in hz)
% outputs:
%     - saccadeDetection - structure containing information about the detected saccades. It also reports the parameters used to make the detections

standardDeviationThreshold = 4; %z threshold for eye velocity
minimumAmplitude = 1.5; %deg
varianceZCuttoff = 4; %zscore for variance cutoff
minimumTimeDelay = 150; %ms

%% Start with a simple thresholding of the eye velocity
eyeVelocity = diff(eyeTrace);
eyeVelocity_normed = eyeVelocity - mean(eyeVelocity); %ensure centering around 0;

candidateSaccades = abs(eyeVelocity_normed) > standardDeviationThreshold*std(eyeVelocity_normed);
saccadeIndices = find(candidateSaccades);

%% Refine the results
%Remove saccades that don't meet a few paramters:
%   1) A steady state plateou must be reached after the saccade that is at
%   least a minimum distance (minimumAmplitude) away from where the saccade
%   was at the beginning of the trace
%   2) There must not be abnormally variance before or after the saccade
%   3) There must not be a missing data point before or after the saccade
%   4) There must be a minimum time interval between 2 adjacent saccades

frameDuration = 1000/frameRate; %ms/frame
plateauPreWindow = 200; %ms
plateauPostWindow = 250; %ms
plateauPreFrames = round(plateauPreWindow/frameDuration); %num frames
plateauPostFrames = round(plateauPostWindow/frameDuration); %num frames
plateauWeights = minimumAmplitude.*linspace(1, 0.2, plateauPostFrames);

%determine the average variance within a window
variancePreWindow = 500; %ms
variancePostWindow = 500; %ms
variancePreFrames = round(variancePreWindow/frameDuration); %num frames
variancePostFrames = round(variancePostWindow/frameDuration); %num frames

varianceDistribution = movvar(eyeVelocity, min(variancePreFrames, variancePostFrames));
meanVariance = mean(varianceDistribution);
stdVariance = std(varianceDistribution);
varianceZ = @(localPoints) abs(meanVariance - var(localPoints))/stdVariance; %z score of local variance

%remove saccades that happen too close to the beginning or end of the trace
tooEarlyFrames = max([plateauPreFrames, variancePreFrames])+2;
tooLateFrames = max([plateauPostFrames, variancePostFrames])+2;
saccadeIndices(saccadeIndices < tooEarlyFrames + 1) = [];
saccadeIndices(saccadeIndices > numel(eyeTrace) - tooLateFrames - 1) = [];

goodSaccades = [];
saccadeAmplitudes = [];
for i = 1:numel(saccadeIndices)
    index_i = saccadeIndices(i);
    
    %1) Check that there are steady state plateaus on either side of the
    %saccade
    preWindowPoints = eyeTrace(index_i-plateauPreFrames-1:index_i-1);
    postWindowPoints = eyeTrace(index_i+1:index_i+plateauPostFrames);
    saccadeAmplitudePerPoint = postWindowPoints - mean(preWindowPoints);
    saccadeDirection = sign(eyeVelocity(index_i)); %+1 for up, -1 for down
    
    %each one of the post points must remain sufficiently far from the mean
    %of the prepoints
    pass = 0;
    for j = 1:numel(postWindowPoints)
        if saccadeDirection*saccadeAmplitudePerPoint(j) < plateauWeights(j)
            pass = 1;
            break
        end
    end
    
    %prepoints shouldn't be too variable:
    if std(preWindowPoints) > 1; pass = 1; end
    
    if pass
        continue
    end
    
    
    %2) check local variance
    preWindowVariancePoints = eyeVelocity(index_i - variancePreFrames-1:index_i-2);
    postWindowVariancePoints = eyeVelocity(index_i + 2:index_i+2+variancePostFrames);
    preVarianceZ = varianceZ(preWindowVariancePoints);
    postVarianceZ = varianceZ(postWindowVariancePoints);
    
    if preVarianceZ > varianceZCuttoff || postVarianceZ > varianceZCuttoff
        continue
    end
    
    %3) check how many whether unknown frames there are across adjacent
    %to the current saccade (velocity should show up as 0 only when there
    %is an unknown)
    percentUnknownsAllowed = 10;
    localVelocities = [preWindowVariancePoints, postWindowVariancePoints];
    if sum(localVelocities == 0) > percentUnknownsAllowed*numel(localVelocities)/100
        continue
    end
    
    %4) Check that this saccade is at least the minimum time interval from
    %the previous
    if numel(goodSaccades) > 0
        lastSaccade = goodSaccades(end);
        numFrames = index_i - lastSaccade;
        interval = numFrames*frameDuration;
        if interval < minimumTimeDelay
            continue
        end
    end
    
    goodSaccades(end+1) = index_i; %if all requirements have been met, add this saccade
end

%% Calculate amplitude, duration, and direction of each detected saccade
goodSaccadeAmplitudes = zeros(1, numel(goodSaccades));
goodSaccadeDirections = zeros(1, numel(goodSaccades));
goodSaccadeDurations = zeros(1, numel(goodSaccades));
maxTimePerSaccade = 200; %ms
maxFramesPerSaccade = round(maxTimePerSaccade/frameDuration);
for i = 1:numel(goodSaccades)
    saccadeDirection = sign(eyeVelocity(goodSaccades(i)));
    
    count = 1;
    while count <= maxFramesPerSaccade
        nextDirection = sign(eyeVelocity(goodSaccades(i)+count));
        if sign(nextDirection) ~= saccadeDirection || count == maxFramesPerSaccade
            saccadeAmplitude = eyeTrace(goodSaccades(i)+count) - eyeTrace(goodSaccades(i));
            break
        else
            count = count + 1;
        end
    end
    
    goodSaccadeAmplitudes(i) = saccadeAmplitude;
    goodSaccadeDurations(i) = count;
    goodSaccadeDirections(i) = saccadeDirection;
end

% %plot saccade locations
% figure
% plot(eyeTrace)
% hold on
% scatter(goodSaccades, eyeTrace(goodSaccades).*ones(1, numel(goodSaccades)), 'dk')
%% Report results in an output structure
%fill in results
saccadeDetection.SaccadeStartFrameNumbers = goodSaccades;
saccadeDetection.SaccadeDurations_frames = goodSaccadeDurations;
saccadeDetection.SaccadeAmplitudes = goodSaccadeAmplitudes;
saccadeDetection.SaccadeDirections = goodSaccadeDirections;

%fill in detection parameters
saccadeDetection.DetectionParameters.velocityZThreshold = standardDeviationThreshold;
saccadeDetection.DetectionParameters.minimumAmplitude = minimumAmplitude;
saccadeDetection.DetectionParameters.varianceZThreshold = varianceZCuttoff;
saccadeDetection.DetectionParameters.minimumTimeDelay = minimumTimeDelay;
saccadeDetection.DetectionParameters.plateauPreWindow = plateauPreWindow;
saccadeDetection.DetectionParameters.plateauPostWindow = plateauPostWindow;
saccadeDetection.DetectionParameters.variancePreWindow = variancePreWindow;
saccadeDetection.DetectionParameters.variancePostWindow = variancePostWindow; 
saccadeDetection.DetectionParameters.maxTimePerSaccade = maxTimePerSaccade;
saccadeDetection.DetectionParameters.addedManually = [];
saccadeDetection.DetectionParameters.removedManually = [];
end