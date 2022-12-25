coverageFactorsL = [2, 3, 6];
numberOfMosaicsL = [1, 2, 3];
for i = 1:numel(coverageFactorsL)
    coverageFactor = coverageFactorsL(i);
    for j = 1:numel(numberOfMosaicsL)
        numberOfMosaics = numberOfMosaicsL(j);
        disp(['Coverage Factor: ' num2str(coverageFactor) ', Mosaics: ' num2str(numberOfMosaics)])
        mosaicSimulation_Master;
    end
end
