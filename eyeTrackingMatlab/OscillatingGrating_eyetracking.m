function [resultsStruct, resultsStruct_batch] = OscillatingGrating_eyetracking(app)
%analysis for the oscillating grating direction stimulus for OKR behavior.

%this function is called by the eyetrackingTracer app (which is the input
%argument)

%can only run after a pupil calibration has been set for this animal.

%Uses data from: .avi recorded video (streampix), .JSON (bassoon export), .CSV
%Pupil data (Deep Lab Cut), and .CSV DAQ (streampix).


%% Run the basics
AnalysisBasics 

%% Get oscillating grating information
numEpochs = stimulusInfo.loggedStimuli.stimulusReps * numel(stimulusInfo.loggedStimuli.gratingOrientations);
if numEpochs ~= numberOfDAQONEpochs || numEpochs ~= numberOfDAQOFFEpochs
    warning('DAQ ANALYSIS DID NOT FIND THE EXPECTED NUMBER OF EPOCHS')
else
    disp('...found the expected number of epochs')
end

preTime = stimulusInfo.loggedStimuli.x_actualPreTime; %s
stimTime = stimulusInfo.loggedStimuli.x_actualStimTime; %s
tailTime = stimulusInfo.loggedStimuli.x_actualTailTime; %s
interStimulusInterval = stimulusInfo.loggedStimuli.x_actualInterStimulusInterval;

%% Get pupil position for every frame
[horizontalPosition, verticalPosition, pupilRadii] = trackPupilPerFrame(pupilData, linearRpMdl, sideDistance_y);

%% Calculate Saccade locations
verticalSaccades = findSaccades(verticalPosition, FrameRate_recording);
horizontalSaccades = findSaccades(horizontalPosition, FrameRate_recording);

%% Divide traces up epoch by epoch 
%calculate the number of frames during each pretime and tailtime
preTimeFrames_record = round(preTime * FrameRate_recording);
tailTimeFrames_record = round(tailTime * FrameRate_recording);

%video frames during the stim time
videoFramesPerEpoch = min(DAQ_OFF - DAQ_ON - preTimeFrames_record - tailTimeFrames_record);

epochStartFrames_absolute = zeros(1, numEpochs); %absolute frame number that this epoch started on
epochEndFrames_absolute = zeros(1, numEpochs);
stimOnFrameByEpoch_absolute = zeros(1, numEpochs); %frame number that marks the start of the stimTime relative - relative to the beginning of this epoch
stimOffFrameByEpoch_absolute = zeros(1, numEpochs);
verticalEyeTraceByEpoch = cell(numEpochs, 1);
horizontalEyeTraceByEpoch = cell(numEpochs, 1);
pupilRadiiByEpoch = cell(numEpochs, 1);

%% determine timing of each epoch based on DAQ data
for i = 1:numEpochs
    verticalEyeTraceByEpoch{i} = verticalPosition(DAQ_ON(i):DAQ_OFF(i));
    horizontalEyeTraceByEpoch{i} = horizontalPosition(DAQ_ON(i):DAQ_OFF(i));
    pupilRadiiByEpoch{i} = pupilRadii(DAQ_ON(i):DAQ_OFF(i));
    
    epochStartFrames_absolute(i) = DAQ_ON(i);
    epochEndFrames_absolute(i) = DAQ_OFF(i);
    
    StartFrame_i = DAQ_ON(i) + preTimeFrames_record; %not taking data associated with pretime and tailtime for now
    EndFrame_i = DAQ_OFF(i) - tailTimeFrames_record;
   
    numFrames_i = EndFrame_i - StartFrame_i;
    if numFrames_i ~= videoFramesPerEpoch
        numExtraFrames = numFrames_i - videoFramesPerEpoch;
        EndFrame_i = EndFrame_i - numExtraFrames; %trim extra frames from the end of the stimulus time
    end
    
    stimOnFrameByEpoch_absolute(i) = StartFrame_i;
    stimOffFrameByEpoch_absolute(i) = EndFrame_i;
end

%% Calculate stimulus velocity on every frame

%Recreate the speed of the grating for each recording frame (assuming
%constant recording frame rate)
cyclesPerStimFrame = stimulusInfo.loggedStimuli.x_numCyclesToShiftByFrame;
secondsPerStimFrame = 1/stimulusInfo.loggedStimuli.x_FR;

%get the time stamp from the onset of the stimulus for every stimulus frame
timeStampsVelocity = secondsPerStimFrame*[0:numel(cyclesPerStimFrame)-1]+secondsPerStimFrame/2;

%Get the time stamp from the start of the stimulus for every recording
%frame
secondsPerRecordFrame = 1/FrameRate_recording;
timeStampsRecord = secondsPerRecordFrame*[0:videoFramesPerEpoch-1]+secondsPerRecordFrame/2;

%match stimulus speeds based on the timing of stimulus and record frames
cyclesPerFrame_recordIndexed = zeros(size(timeStampsRecord));
for i = 1:videoFramesPerEpoch
    currentRecordTime = timeStampsRecord(i);
    [~, nearestStimTime] = min(abs(timeStampsVelocity - currentRecordTime));
    cyclesPerFrame_recordIndexed(i) = cyclesPerStimFrame(nearestStimTime);
end

%transform cycles per record frame to visual degrees per second
gratingSpeedPerFrame = cyclesPerFrame_recordIndexed * stimulusInfo.loggedStimuli.x_FR...
    / stimulusInfo.loggedStimuli.spatialFrequency; % this is speed of the grating on each recording frame in deg/s

recordFramesPerOscillation = stimulusInfo.loggedStimuli.oscillationPeriod *FrameRate_recording;

%determine number of frames after the start of the pretime that the
%oscillation gets back to the beginning of its cycle (i.e. first time
%stimulus speed = 0 and oscillation is moving in the correct direction)
firstOscillationStartOffset = (recordFramesPerOscillation/4)*(stimulusInfo.loggedStimuli.oscillationPhaseShift/90);


%% Save to the final structures
resultsStruct.StimFile = stimulusF;
resultsStruct.PupilFile = pupilF;
resultsStruct.DAQFile = DAQF;
resultsStruct.stimulusInfo = stimulusInfo;
resultsStruct.meta = stimulusInfo.loggedStimuli;
resultsStruct.protocolName = protocolName; %saves the name of the protocol alongside the calibration
resultsStruct.DAQ_ON = DAQ_ON;
resultsStruct.DAQ_OFF = DAQ_OFF;
resultsStruct.FrameRate_recording = FrameRate_recording;
resultsStruct.FramesPerEpoch_recording = videoFramesPerEpoch;
resultsStruct.AllOrientations = stimulusInfo.loggedStimuli.gratingOrientations;
resultsStruct.orientationByEpoch = stimulusInfo.loggedStimuli.x_orientationLog;
resultsStruct.numEpochs = numEpochs;
resultsStruct.preTime = preTime;
resultsStruct.stimTime = stimTime;
resultsStruct.tailTime = tailTime;
resultsStruct.interStimulusInterval = interStimulusInterval;
resultsStruct.FrameHeight_recording = frameHeight;
resultsStruct.FrameWidth_recording = frameWidth;
resultsStruct.epochStartFrames_absolute = epochStartFrames_absolute; %absolute frame number that this epoch started on
resultsStruct.epochEndFrames_absolute = epochEndFrames_absolute;
resultsStruct.stimOnFrameByEpoch_absolute = stimOnFrameByEpoch_absolute; %relative frame number (from the start of the epoch) that the stimulus turned on at
resultsStruct.stimOffFrameByEpoch_absolute = stimOffFrameByEpoch_absolute;
resultsStruct.horizontalSaccades = horizontalSaccades;
resultsStruct.verticalSaccades = verticalSaccades;
resultsStruct.recordingFramesPerOscillationPeriod = recordFramesPerOscillation;
resultsStruct.OscillationPhaseShift_frames = firstOscillationStartOffset;

resultsStruct_batch = resultsStruct; %copy over data to resultsStruct_batch before including full eye traces

resultsStruct.verticalEyeTraceByEpoch = verticalEyeTraceByEpoch;
resultsStruct.horizontalEyeTraceByEpoch = horizontalEyeTraceByEpoch;
resultsStruct.pupilRadiiByEpoch = pupilRadiiByEpoch;
resultsStruct.verticalEyeTrace_full = verticalPosition;
resultsStruct.horizontalEyeTrace_full = horizontalPosition;
resultsStruct.pupilRadii_full = pupilRadii;
resultsStruct.gratingSpeedPerFrame = gratingSpeedPerFrame; %in deg/s the speed of the stimulus on each recording frame of the stim time

%mark down that the fields that exist in the full structure but are only
%referenced in the batch structure (i.e. these are the data-heavy fields)
resultsStruct_batch.verticalEyeTraceByEpoch = 'exists';
resultsStruct_batch.horizontalEyeTraceByEpoch = 'exists''exists';
resultsStruct_batch.pupilRadiiByEpoch = 'exists';
resultsStruct_batch.verticalEyeTrace_full = 'exists';
resultsStruct_batch.horizontalEyeTrace_full = 'exists';
resultsStruct_batch.pupilRadii_full = 'exists';
resultsStruct_batch.gratingSpeedPerFrame = 'exists';

%% Make the figure
verticalSaccadeTimePoints = verticalSaccades.SaccadeStartFrameNumbers;
horizontalSaccadeTimePoints = horizontalSaccades.SaccadeStartFrameNumbers;

%concatonate vectors to recreate the stimulus speed on each frame (0 during
%pre, post and isi time
fullStimulusRecreation = zeros(1, numel(verticalPosition));
for i = 1:numel(epochStartFrames_absolute)
    fullStimulusRecreation(stimOnFrameByEpoch_absolute(i):...
        stimOffFrameByEpoch_absolute(i) - 1) = gratingSpeedPerFrame;
end


figure
title([app.animalName,': ', protocolName])
hold on
plot(verticalPosition,'Color', '#0072BD')
s1 = scatter(verticalSaccadeTimePoints, verticalPosition(verticalSaccadeTimePoints), 'kd');
s1.HandleVisibility = 'off';
plot(horizontalPosition,'Color', '#D95319')
s2 = scatter(horizontalSaccadeTimePoints, horizontalPosition(horizontalSaccadeTimePoints), 'kd');
s2.HandleVisibility = 'off';
plot(pupilRadii, 'Color', '#77AC30')
plot(fullStimulusRecreation, 'Color', '#7E2F8E')
allDAQs = [DAQ_ON, DAQ_OFF];
for i = 1:numel(allDAQs)
    p = plot([allDAQs(i) allDAQs(i)], [-30 50], '--r');
    p.HandleVisibility = 'off';
end
legend('Vertical Pos.', 'Horizontal Pos.', 'Pupil Radius', 'Stimulus')
ylabel('Degrees')
xlabel('Frame Number')
end