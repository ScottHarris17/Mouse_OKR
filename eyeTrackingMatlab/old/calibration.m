[pupil_fname, path] = uigetfile('.csv', 'Select Pupil Data');

pupilLabelIndex = strfind(pupil_fname, '_Pupil.');

DAQ_fname = [pupil_fname(1:pupilLabelIndex) 'DAQ.csv'];

pupilFile = fullfile(path, pupil_fname);
pupilData = readmatrix(pupilFile);

DAQfile = fullfile(path, DAQ_fname);
DAQData = readmatrix(DAQfile);

slashes = strfind(path, '\'); %change for mac/unix
expName = path(slashes(end-1)+1:slashes(end)-1);
animalName_temp = path(slashes(end-2)+1:slashes(end-1)-1);
animalName = animalName_temp(1:strfind(animalName_temp, '_record')-1);
JSONPath = [path(1:slashes(end-2)), animalName, '_stimulus\', expName '_Complete.json'];

% [JSONf, JSONp] = uigetfile('.JSON', 'Select JSON Data');
% JSONPath = fullfile(JSONp, JSONf);

%Get JSON data
fid = fopen(JSONPath); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
stimulusInfo = jsondecode(str);


%% Basic information
frameHeight = 296; %pix
frameWidth = 512; %pix

numFrames = size(pupilData, 1);

%flip y coordinates because image coordinate system is upside down
pupilData(:, [3:3:27]) = frameHeight - pupilData(:, [3:3:27]);

%DAQ ON frames
DAQ_diff = diff(DAQData(:, 2));
DAQ_ON = find(DAQ_diff > 2.5); %when the DAQ command switches on
DAQ_OFF = find(DAQ_diff < -2.5); %when the DAQ command switches off

%light level sequence
lightLevels_all = stimulusInfo.loggedStimuli.x_lightLevelLog;

if isfield(stimulusInfo.loggedStimuli, 'x_versionNumber')
else
    %for version 1.1 of the PupilCalibration Bassoon stimulus, the
    %lightLevelLog does not accurately reflect the order of stimuli that
    %were actually played. Instead, the first stimulus was always full dark
    %(value of -1), then the stimuli played as listed in the lightLevelLog
    %up to the very last stimulus, which wasn't played. This version of the
    %stimulus also didn't have an explicit versionNumber attribute saved
    %with it, so animals for which this stimulus was used wouldn't have
    %this field in the log.
    lightLevels_all = [-1; lightLevels_all(1:end-1)]; %version 1.1 of the PupilCalibration Basso
end


%angle swing of camera during calibration
angleSwing = stimulusInfo.loggedStimuli.angle; %in deg
angleSwing_rad = deg2rad(angleSwing);

%recording date
recordingDate = datestr(stimulusInfo.experimentStartTime);

%% Load the video
%vidObj = VideoReader('C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\behavior\Eyetracking\SHOKRB1\SHOKRB1_record\PupilCalibration20220113\SHOKRB1_PupilCalibration20220113.avi');
%% Find pupil size and location for each DAQ frame

distanceAndPupilSize = zeros(numel(lightLevels_all), 2);
for i = 1:2:numel(lightLevels_all)*2    
    ONFrameLeft = DAQ_ON(i); %camera should swing left first
    OFFFrameLeft = DAQ_OFF(i);
    
    ONFrameRight = DAQ_ON(i+1); %camera should swing right second
    OFFFrameRight = DAQ_OFF(i+1);
    
%     videoFrames = squeeze(read(vidObj, [ONFrameLeft, OFFFrameLeft]));
%     imshow(videoFrames(:, :, 1));
    
    frameNums_Left = ONFrameLeft:OFFFrameLeft;
    frameNums_Right = ONFrameRight:OFFFrameRight;
    
    xDistance_Left = []; %keeps track of difference between pupil center and CR
    radii_Left = [];
    for j = 1:numel(frameNums_Left)
        trackingData_j = pupilData(frameNums_Left(j), :);
        [pPosition, pRadius, CR] = findPupilPosition(trackingData_j, 0.5);
        if pPosition == -1
            continue
        end
        
        xDistance_Left(end+1) = pPosition(1) - CR(1); %Calculate horizontal distance between pupil and corneal reflection
        radii_Left(end+1) = pRadius;
    end
    
    xDistance_Right = []; %keeps track of difference between pupil center and CR
    radii_Right = [];
    for j = 1:numel(frameNums_Right)
        trackingData_j = pupilData(frameNums_Right(j), :);
        [pPosition, pRadius, CR] = findPupilPosition(trackingData_j, 0.5);
        if pPosition == -1 %bad frame
            continue
        end
        
        xDistance_Right(end+1) = pPosition(1) - CR(1); %Calculate horizontal distance between pupil and corneal reflection
        radii_Right(end+1) = pRadius;
    end
    
    if numel(xDistance_Left) <= 2 || numel(xDistance_Right) <= 2
        disp('No good frames')
        distanceAndPupilSize((i+1)/2, :) = [NaN, NaN];
        continue
    end

    Rp = abs((mean(xDistance_Left)-mean(xDistance_Right))/angleSwing_rad);
    pupilSize = mean([radii_Left, radii_Right]);
    
    distanceAndPupilSize((i+1)/2, :) = [Rp, pupilSize];
end

lightLevels_unique = unique(lightLevels_all);
meanDistanceByLightLevel = zeros(1, numel(lightLevels_unique));
meanPupilSizeByLightLevel = zeros(1, numel(lightLevels_unique));
for i = 1:numel(lightLevels_unique)
    indices = lightLevels_all == lightLevels_unique(i);
    
    distances = distanceAndPupilSize(indices, 1);
    pupilSizes = distanceAndPupilSize(indices, 2);
    
    distances_valid = distances(~isnan(distances));
    pupilSizes_valid = pupilSizes(~isnan(pupilSizes));
    
    if numel(distances_valid) == 0
        warning(['No valid calibrations for light level ' num2str(lightLevels_unique(i))])
        meanDistanceByLightLevel(i) = NaN;
        meanPupilSizeByLightLevel(i) = NaN;
        continue
    end
    
    meanDistanceByLightLevel(i) = mean(distances_valid);
    meanPupilSizeByLightLevel(i) = mean(pupilSizes_valid);
end

%Generate the linear model for predicting Rp from pupil size
[R, p_value] = corrcoef(meanPupilSizeByLightLevel, meanDistanceByLightLevel);
Slope = R(2) * std(meanDistanceByLightLevel)/std(meanPupilSizeByLightLevel);
Intercept = mean(meanDistanceByLightLevel) - Slope*mean(meanPupilSizeByLightLevel);
linearRpMdl = @(pupilSize) Slope*pupilSize + Intercept; %function for predicting Rp from pupil size

figure
title('Pupil size and Distance')
hold on
scatter(meanPupilSizeByLightLevel, meanDistanceByLightLevel, 'k', 'filled')
plot(0:100, linearRpMdl(0:100), '--r')
xlabel('Pupil Radius (px)')
ylabel('Distance')
legend({'Observations', 'Model'})

%% Determine location of pupil equator with side LED
topLED_ON = DAQ_ON(end-1);
topLED_OFF = DAQ_OFF(end -1);

sideLED_ON = DAQ_ON(end);
sideLED_OFF = DAQ_OFF(end);

topLED_frames = topLED_ON:topLED_OFF;
sideLED_frames = sideLED_ON:sideLED_OFF;

topCRs_x = [];
topCRs_y = [];
for i = 1:numel(topLED_frames)
    trackingData = pupilData(topLED_frames(j), :);
    [pPosition, pRadius, CR] = findPupilPosition(trackingData, 0.5);
    if pPosition == -1 %bad frame
        continue
    end
    topCRs_x(end + 1) = CR(1);
    topCRs_y(end + 1) = CR(2);
end
if numel(topCRs_x) == 0
    warning('Could not resolve corneal reflection for top led')
end


sideCRs_x = [];
sideCRs_y = [];
for i = 1:numel(topLED_frames)
    trackingData = pupilData(sideLED_frames(j), :);
    [pPosition, pRadius, CR] = findPupilPosition(trackingData, 0.5);
    if pPosition == -1 %bad frame
        continue
    end
    sideCRs_x(end + 1) = CR(1);
    sideCRs_y(end + 1) = CR(2);
end
if numel(sideCRs_x) == 0
    warning('Could not resolve corneal reflection for side led')
end

sideDistance_x = mean(sideCRs_x) - mean(topCRs_x);
sideDistance_y = mean(sideCRs_y) - mean(topCRs_y);
 clearvars('-except', 'linearRpMdl', 'sideDistance_x', 'sideDistance_y')