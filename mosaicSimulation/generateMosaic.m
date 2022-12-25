function centers = generateMosaic(coverageFactor, cellRadius, noise, retinaRadius, varargin)
% This function generates model RGC mosaics that follow the following
% rules: 1) no two cells can overlap, regardless of if they are members of the
% same mosaic, 2) cells that are members of the same mosaic obey a (noisey)
% exclusion zone. Flat mount retinas are modeled as a circle. RGC density
% is modeled as an exponential function of eccentricity.
% Inputs:
%     - coverageFactor: Average cells per unit area
%     - cellRadius: radius of the dendritic tree of each cell
%     - noise: noise scale (see below for implementation)
%     - retinaRadius: radius of the model retina
%     - varargin: nx2 matrix holding x and y coordinates of cells from other mosaics
% Outputs:
%     - centers: nx2 matrix holding the x and y coordinates of each cell in a complete mosaic

numArgs = numel(varargin);
if numArgs > 1
    otherCenters = varargin{2};
else
    otherCenters = [];
end
otherCentersMinDistance = 15;%um, for cells that are part of other mosaics

%Calculated Retina Parameters
retinaArea = pi*retinaRadius^2;
cellArea = pi*cellRadius^2;
numCells = ceil(coverageFactor*retinaArea/cellArea);

%% Initialize the retina and RGCs
centers = [];
minDistance = cellRadius/(sqrt(coverageFactor)); %exclusion zone as a function of coverage factor and cell size
for i = 1:numCells
    checkResult = 0;
    noiseAdjustment = normrnd(0, noise)*minDistance; %some gaussian jitter is added.
    while ~checkResult
        [theta, rho] = pickCoordinates(retinaRadius); %see function below
        [x, y] = pol2cart(theta, rho);
        checkResult = checkDistances(x, y, centers, minDistance, noiseAdjustment, otherCenters, otherCentersMinDistance); %see function below
    end
    centers(i, 1) = x;
    centers(i, 2) = y;
end

showFigure = 1;
if nargin > 3
    for i = 1:numel(varargin)
        arg_i = varargin{i};
        if ischar(arg_i) && strcmp(arg_i, 'SupressFigure')
            showFigure = 0;
        end
    end
end

if showFigure
    figure()
    scatter(centers(:, 1), centers(:, 2), 'filled', 'MarkerFaceColor', [175, 107, 51]./255,  'MarkerFaceAlpha', 0.41)
end


function [theta, rho] = pickCoordinates(radius)
    %generates random polar coordinates normalized to radius size. Using
    %random polar coordinates will scale density by eccentricity to roughly
    %approximate the physiological dependency of RGC densities on
    %eccentricity
    theta = rand()*2*pi; %random angle in radians
    rho = rand()*radius; %random radius
end

function checkResult = checkDistances(x, y, cellLocations, minDistance, noiseAdjustment, otherCenters, otherMinDistance)
%     checks that the location of a new cell does not violate exclusion zone
%     of another cell based on the minumum distance specified.
%     Inputs:
%         - x: x coordinate of the cell
%         - y: y coordinate of the cell
%         - cellLocations: nx2 matrix holding x,y coordiantes of all other cells that are members of the same mosaic
%         - minDistance: minimum distance that can exist between two cells of the same mosaic
%         - noiseAdjustment: random noise for the current cell
%         - otherCenters: nx2 matrix holding x,y coordinates of cells from other mosaics
%         - otherMinDistance: minimum distance that can exist between two cells of different mosaics.
%      Outputs:
%         -checkResult: bool. 1 if the new cell's coordinates obey the rules. Otherwise 0.
    checkResult = 1; 
    for n = 1:size(cellLocations, 1) %first check for the same mosaic
        existantX = cellLocations(n, 1);
        existantY = cellLocations(n, 2);
        distance = sqrt((x - existantX)^2 + (y-existantY)^2);
        if distance < minDistance + noiseAdjustment %some jitter added
            checkResult = 0;
            return
        end           
    end
    
    %Can't have overlapping cells regardless of type... if cells of another
    %type exist within this retina, check that new cell does not overlap
    %with it
    if numel(otherCenters)
        for m = 1:size(otherCenters, 1)
                otherX = otherCenters(m, 1);
                otherY = otherCenters(m, 2);
                distance = sqrt((x - otherX)^2 + (y-otherY)^2);
                if distance < otherMinDistance
                    checkResult = 0;
                    return
                end
        end
    end
end
end
