%% This script is called by eyetracking analysis functions. It runs the basic things
% that these functions have in common - such as loading the data files and
% determining the basic information about them

%% load data
selection = app.SelectionIndexEditField.Value;

protocolName = app.UITable.Data{selection, 1};
stimulusF = app.UITable.Data{selection, 2};
pupilF = app.UITable.Data{selection, 3};
DAQF = app.UITable.Data{selection, 4};

pupilData = readmatrix(pupilF);
DAQData = readmatrix(DAQF);

%Get JSON data
fid = fopen(stimulusF); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
stimulusInfo = jsondecode(str);

%% Get calibration data
linearRpMdl = app.eyetrackingParameters.Calibrations.(app.activeCalibration).linearModel.mdl;
sideDistance_y = app.eyetrackingParameters.Calibrations.(app.activeCalibration).sideLEDDistance_y;
resultsStruct.CalibrationInfo.CalibrationName = app.activeCalibration;
resultsStruct.CalibrationInfo.linearRpMdl = linearRpMdl;
resultsStruct.CalibrationInfo.sideDistance_y = sideDistance_y;

%% Basic information
%camera info
frameHeight = 296; %pix
frameWidth = 512; %pix

numFrames = size(pupilData, 1);

%flip y coordinates because image coordinate system is upside down
pupilData(:, [3:3:27]) = frameHeight - pupilData(:, [3:3:27]);

%DAQ ON frames
DAQ_diff = diff(DAQData(:, 2));
allDAQChanges = find(abs(DAQ_diff) > 2.5);
lastDAQ = 0;
adjustedDAQData = DAQData;
minimumDAQInterval = 4; %minimum number of frames between DAQ switches
for i = 1:numel(allDAQChanges)
    thisDAQ = allDAQChanges(i);
    if thisDAQ - lastDAQ < minimumDAQInterval
        adjustedDAQData(lastDAQ + 1:thisDAQ, 2) = DAQData(lastDAQ - 1, 2);
    end
    lastDAQ = thisDAQ;
end

adjustedDAQ_diff = diff(adjustedDAQData(:, 2));
DAQ_ON = find(adjustedDAQ_diff > 2.5); %when the DAQ command switches on
DAQ_OFF = find(adjustedDAQ_diff < -2.5); %when the DAQ command switches off

numberOfDAQONEpochs = numel(DAQ_ON);
numberOfDAQOFFEpochs = numel(DAQ_OFF);

VideoInfoFile = app.UITable.Data{selection, 5};
if contains(VideoInfoFile, '.avi')
    disp('Downloading video file')
    VidObj = VideoReader(app.UITable.Data{selection, 5});
    FrameRate_recording = VidObj.FrameRate;
    disp('Saving frame rate to .mat file for future access')
    vidInfo.FrameRate = FrameRate_recording;
    save(fullfile(app.animalPath, app.UITable.Data{selection, 1}, 'videoInfo.mat'), 'vidInfo')
else
    load(VideoInfoFile, 'vidInfo'); %loads vidInfo.mat for this stimulus
    FrameRate_recording = vidInfo.FrameRate;
end