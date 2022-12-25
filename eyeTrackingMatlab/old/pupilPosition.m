path = 'C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\behavior\Eyetracking\SHOKRB1\SHOKRB1_record';
fname = '\MovingGratingUpAndDown\SHOKRB1_MovingGratingDirection_UpAndDown_Pupil';
f = fullfile(path, fname);
m = readmatrix(f);

%Get JSON data
[fname, path] = uigetfile('.json'); %load '_completed_modified'
f = fullfile(path, fname);
fid = fopen(f); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
stimulusInfo = jsondecode(str);

%%
frameHeight = 296; %pix
frameWidth = 512; %pix

numFrames = size(m, 1);

%flip y coordinates because image coordinate system is upside down
m(:, [3:3:27]) = frameHeight - m(:, [3:3:27]);


LikelihoodThreshold = 0.5;
lastKnownVals = NaN;
knownFrames = ones(1, numFrames);
xDisplacement = zeros(1, numFrames);
yDisplacement = zeros(1, numFrames);
pupilRadius = zeros(1, numFrames);
g = waitbar(0, 'Please Wait');
for i = 1:numFrames
    waitbar(i/numFrames, g, 'Please Wait')
    thisFrame = m(i, :);
    thisFrameVals = struct();
    thisFrameVals.pN = thisFrame([2, 3]);
    thisFrameVals.pNE = thisFrame([5, 6]);
    thisFrameVals.pE = thisFrame([8,9]);
    thisFrameVals.pSE = thisFrame([11, 12]);
    thisFrameVals.pS = thisFrame([14, 15]);
    thisFrameVals.pSW = thisFrame([17, 18]);
    thisFrameVals.pW = thisFrame([20, 21]);
    thisFrameVals.pNW = thisFrame([23, 24]);
    thisFrameVals.CR = thisFrame([26, 27]);
    thisFrameVals.allLikelihoods = thisFrame(4:3:28);
   
    if min(thisFrameVals.allLikelihoods) < LikelihoodThreshold
        knownFrames(i) = 0; %mark down that you didn't know this frame
        if ~isstruct(lastKnownVals)
            continue %on the event that the first frame has unknown vals, wait until the first valid frame shows up
        end
        thisFrameVals = lastKnownVals;
    else
        lastKnownVals = thisFrameVals;
    end
    
    %calculate pupil center
    linearEqs = {};
    [c1, s1] = calculateLine(thisFrameVals.pN, thisFrameVals.pS);
    linearEqs(1, [1, 2]) = [{c1}, {s1}];
    [c2, s2] = calculateLine(thisFrameVals.pNE, thisFrameVals.pSW);
    linearEqs(2, [1, 2]) = [{c2}, {s2}];
    [c3, s3] = calculateLine(thisFrameVals.pE, thisFrameVals.pW);
    linearEqs(3, [1, 2]) = [{c3}, {s3}];
    [c4, s4] = calculateLine(thisFrameVals.pSE, thisFrameVals.pNW);
    linearEqs(4, [1, 2]) = [{c4}, {s4}];
    
    combinations = [1, 2; 1, 3; 1, 4; 2, 3; 2, 4; 3, 4];
    numCombinations = size(combinations, 1);
    allCenterGuesses = zeros(numCombinations, 2);
    for j = 1:numCombinations
        indices = combinations(j, :);
        A1 = linearEqs{indices(1), 1};
        B1 = linearEqs{indices(1), 2};
        A2 = linearEqs{indices(2), 1};
        B2 = linearEqs{indices(2), 2};
        intersection = linsolve([A1;A2], [B1;B2]);
        allCenterGuesses(j, :) = intersection';
    end
    
    pupilCenter = fliplr(mean(allCenterGuesses));
    
    xDisplacement(i) = pupilCenter(1) - thisFrameVals.CR(1);
    yDisplacement(i) = pupilCenter(2) - thisFrameVals.CR(2);
    
    %calculate pupil size
    radii = zeros(1, 8);
    radii(1) = cartDist(pupilCenter, thisFrameVals.pN);
    radii(2) = cartDist(pupilCenter, thisFrameVals.pNE);
    radii(3) = cartDist(pupilCenter, thisFrameVals.pE);
    radii(4) = cartDist(pupilCenter, thisFrameVals.pSE);
    radii(5) = cartDist(pupilCenter, thisFrameVals.pS);
    radii(6) = cartDist(pupilCenter, thisFrameVals.pSW);
    radii(7) = cartDist(pupilCenter, thisFrameVals.pW);
    radii(8) = cartDist(pupilCenter, thisFrameVals.pNW);
    
    pupilRadius(i) = mean(radii);
end
close(g)
plot(yDisplacement)
        
% scatter(pN(1), pN(2))
% hold on
% scatter(pNE(1), pNE(2))
% scatter(pE(1), pE(2))
% scatter(pSE(1), pSE(2))
% scatter(pS(1), pS(2))
% scatter(pSW(1), pSW(2))
% scatter(pW(1), pW(2))
% scatter(pNW(1), pNW(2))
% scatter(CR(1), CR(2), 'filled');
% scatter(pupilCenter(2), pupilCenter(1), 'filled')
