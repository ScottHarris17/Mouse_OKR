function [pPosition, pRadius, CR] = findPupilPosition(DLC_data, likelihoodThreshold)
%Given deeplabcut output in vector form, this function returns pupil center
%position, pupil radius, and corneal reflection location in x, y
%coordinates on a single frame
%inputs:
%   - DLC_data = single row of the deeplabcut table (/.csv) output as a row
%   vector
%   - likelihoodThreshold = threshold at which to discard data if feature
%   likelihood is less than

%based off of deeplabcut model with 8 features around the pupil and one
%corneal reflection feature
pN = DLC_data([2, 3]);
pNE = DLC_data([5, 6]);
pE = DLC_data([8,9]);
pSE = DLC_data([11, 12]);
pS = DLC_data([14, 15]);
pSW = DLC_data([17, 18]);
pW = DLC_data([20, 21]);
pNW = DLC_data([23, 24]);
CR = DLC_data([26, 27]);
allLikelihoods = DLC_data(4:3:28);

if min(allLikelihoods) < likelihoodThreshold
    pPosition = -1;
    pRadius = -1;
    CR = -1;
    return
end


%calculate pupil center
linearEqs = {};
[c1, s1] = calculateLine(pN, pS);
linearEqs(1, [1, 2]) = [{c1}, {s1}];
[c2, s2] = calculateLine(pNE, pSW);
linearEqs(2, [1, 2]) = [{c2}, {s2}];
[c3, s3] = calculateLine(pE, pW);
linearEqs(3, [1, 2]) = [{c3}, {s3}];
[c4, s4] = calculateLine(pSE, pNW);
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

pPosition = fliplr(mean(allCenterGuesses)); %consider changing to medioid
    
%calculate pupil radius
radii = zeros(1, 8);
radii(1, 1) = cartDist(pPosition, pN);
radii(1, 2) = cartDist(pPosition, pNE);
radii(1, 3) = cartDist(pPosition, pE);
radii(1, 4) = cartDist(pPosition, pSE);
radii(1, 5) = cartDist(pPosition, pS);
radii(1, 6) = cartDist(pPosition, pSW);
radii(1, 7) = cartDist(pPosition, pW);
radii(1, 8) = cartDist(pPosition, pNW);

pRadius = mean(radii);


end

% 
% scatter(pN(1), pN(2))
% hold on
% scatter(pNE(1), pNE(2))
% scatter(pE(1), pE(2))
% scatter(pSE(1), pSE(2))
% scatter(pS(1), pS(2))
% scatter(pSW(1), pSW(2))
% scatter(pW(1), pW(2))
% scatter(pNW(1), pNW(2))
% 
% scatter(CR(1), CR(2), 'filled')
% scatter(pPosition(1), pPosition(2), 'filled')
