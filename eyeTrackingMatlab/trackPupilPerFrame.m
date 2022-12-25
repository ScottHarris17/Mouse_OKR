function [horizontalPosition, verticalPosition, pupilRadii] =...
    trackPupilPerFrame(pupilData, linearRpMdl, sideDistance_y)
%This function calculates the position of the pupil on every frame of a
%video. Requries calibration data.
%Inputs:
%     -pupilData = deeplabcut output (originally in the form of a csv file, converted to matlab table using tableread())
%     -linearRpMdl = linear model for Rp size from calibration
%     -sideDistance_y = distance from the top LED to the side LED in the y direction
% 
%Outputs: 
%     -horizontalPosition = x coordinate of pupil across all frames (in degrees)
%     -verticalPosition = y coordinate of pupil across all frames (in degrees) 
%     -pupilRadii = pupil radius size across frames (in degrees)
CertaintyThreshold = 0.9;

numFrames = size(pupilData, 1);
horizontalPosition = zeros(1, numFrames);
verticalPosition = zeros(1, numFrames);
pupilRadii = zeros(1, numFrames);

lastValidDetection = [123.456, 123.456, 123.456, 123.456, 123.456];
validDetections = ones(1, numFrames);
for i = 1:numFrames
    trackingData = pupilData(i, :);
    [pPosition, pRadius, CR] = findPupilPosition(trackingData, CertaintyThreshold);
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
    pupilRadii(i) = pRadius;
end

end
