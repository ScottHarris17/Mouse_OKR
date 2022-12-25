coverageFactor = 3; %average cells/micron^2 2, 3, 6
dendriticTreeRadius = 141; %Microns
sigma = 0.2; %Noise value (unitless)
retinaRadius = 2400; %um
numberOfMosaics = 1; %integer value

thresholdDistances = 0:10:400; %um - for fraction of pairs that are triples analysis
numDRPBins = 100; %for DRP analysis

%single mosaic info:
probe = generateMosaic(coverageFactor, dendriticTreeRadius, sigma, retinaRadius); %, 'SupressFigure'); %test mosaic generator to understand how many cells will be used
cellsPerMosaic = size(probe, 1);


%simulation repeats
numRuns = 30;

nnDistances = zeros(numRuns, cellsPerMosaic*numberOfMosaics);
nnDistances3 = zeros(size(nnDistances));
drps = zeros(numRuns, numDRPBins); %second dimension is number of bins
fractionTriples = zeros(numRuns, numel(thresholdDistances));
ppm = ParforProgressbar(numRuns);
parfor i = 1:numRuns
    ppm.increment();
    
    %generate the simulated retina
    centers = generateMosaic(coverageFactor, dendriticTreeRadius, sigma, retinaRadius, 'SupressFigure');
    for j = 1:numberOfMosaics - 1
        centers = [centers; generateMosaic(coverageFactor, dendriticTreeRadius, sigma, retinaRadius, 'SupressFigure', centers)];
    end
    
    %do the analyses
    nnDistances(i, :) = nearestNeighbor_SIMULATION(centers);
    [fractionTriples(i, :), nnDistances3(i, :)] = threewayNearestNeighbor_SIMULATION(centers, thresholdDistances);
    drps(i, :) = drp_SIMULATION(centers, numDRPBins, retinaRadius);
end
delete(ppm)

nnDistances(nnDistances == 0)  = nan;
nnDistances3(nnDistances3 == 0) = nan;

nnFig = figure;
histogram(nnDistances(:), 'Normalization', 'probability')
hold on
title(['Nearest Neighbor Distances; Mosaics = ' num2str(numberOfMosaics)...
    '; CPM = ' num2str(cellsPerMosaic)])
xlabel('um')
ylabel('Fraction')
xlim([0 max(thresholdDistances)])


nn3Fig = figure;
histogram(nnDistances3(:), thresholdDistances, 'Normalization', 'probability')
hold on
title(['Nearest Threeway Neighbor Distances; Mosaics = ' num2str(numberOfMosaics)...
    '; CPM = ' num2str(cellsPerMosaic)])
xlabel('Distance um')
ylabel('Fraction')
xlim([0 max(thresholdDistances)])

m3F = mean(fractionTriples);
e3F = std(fractionTriples);
nn3FractionFig = figure;
hold on
% for i = 1:numel(thresholdDistances)
%     scatter(repmat(thresholdDistances(i), numel(fractionTriples(:, i)), 1), fractionTriples(:, i), 'MarkerEdgeColor', [0.2, 0.2, 0.2])
% end
shadedErrorBar(thresholdDistances, m3F, e3F)
title(['Fraction of Pairs that are In Triples; Mosaics = ' num2str(numberOfMosaics)...
    '; CPM = ' num2str(cellsPerMosaic)])
xlabel('Distance From Reference Cells (um)')
ylabel('Fraction of Pairs')
ylim([0 1])
hold off

drpFig = figure;
title(['Density Recovery Profile; Mosaics = ' num2str(numberOfMosaics)...
    '; CPM = ' num2str(cellsPerMosaic)])
hold on
box off
semilogx(logspace(1, log10(retinaRadius*2), numDRPBins)./1000, mean(drps)); %in pixels
set(gca, 'XScale', 'log');
xlabel('mm')
ylabel('cells/mm^2')
plot(logspace(1, log10(retinaRadius*2), numDRPBins)./1000,...
    repmat(cellsPerMosaic*numberOfMosaics/(pi*(retinaRadius/1000)^2), 1, numDRPBins), '--r');


ROOT = 'C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\Images\MTN\mosaicSimulations';

saveas(nnFig, [ROOT, '/figures/sigma2/mosaics_' num2str(numberOfMosaics) '_coverage_'...
    num2str(coverageFactor) '_nearestNeighbor.svg'])
saveas(nn3Fig, [ROOT '/figures/sigma2/mosaics_' num2str(numberOfMosaics) '_coverage_'...
    num2str(coverageFactor) '_nearestNeighbor3.svg'])
saveas(nn3FractionFig, [ROOT '/figures/sigma2/mosaics_' num2str(numberOfMosaics) '_coverage_'...
    num2str(coverageFactor) '_fractionTriples.svg'])
saveas(drpFig, [ROOT '/figures/sigma2/mosaics_' num2str(numberOfMosaics) '_coverage_'...
    num2str(coverageFactor) '_DRP.svg'])
close all