%% run the pupil calibration
calibration

%% Load data for actual stimulus
[pupil_fname, path] = uigetfile('.csv', 'Find pupil position');

pupilFile = fullfile(path, pupil_fname);
pupilData = readmatrix(pupilFile);

frameHeight = 296; %pix
frameWidth = 512; %pix

numFrames = size(pupilData, 1);

%flip y coordinates because image coordinate system is upside down
pupilData(:, [3:3:27]) = frameHeight - pupilData(:, [3:3:27]);

%%
horizontalPosition = zeros(1, numFrames);
verticalPosition = zeros(1, numFrames);

lastValidDetection = [123.456, 123.456, 123.456, 123.456, 123.456];
validDetections = ones(1, numFrames);
g = waitbar(0, 'Please Wait');
for i = 1:numFrames
    waitbar(i/numFrames, g, 'Please Wait')
    trackingData = pupilData(i, :);
    [pPosition, pRadius, CR] = findPupilPosition(trackingData, 0.5);
    if pPosition == -1
        pPosition = lastValidDetection(1:2); pRadius = lastValidDetection(3); CR = lastValidDetection(4:5);
        validDetections(i) = 0; %mark that this was not a validframe
    else
        lastValidDetection = [pPosition, pRadius, CR];
    end
    
    dX = pPosition(1) - CR(1);
    dY_Rp = pPosition(2) - CR(2);
    dY = dY_Rp - sideDistance_y;
    
    Rp = linearRpMdl(pRadius);
    
    Rp0 = sqrt(Rp.^2 + dY.^2);
    horizontalAngle = rad2deg(asin(dX/Rp));
    verticalAngle = rad2deg(asin(dY/Rp0));
    
    horizontalPosition(i) = horizontalAngle;
    verticalPosition(i) = verticalAngle;
end
close(g)
% plot(verticalPosition)
