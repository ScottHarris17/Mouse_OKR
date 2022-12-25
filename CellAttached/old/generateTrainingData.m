function trainingData = generateTrainingData(stimuli, runsPerStimulus, cdfs, varargin)
%inputs:
%stimuli - 1d vector of the stimuli orientations (or categories) to use
%runsPerStimulus - integer value, number of time to simulate each stimulus
%cdfs - 1xn cell array of CDFs for distributions of cells to get responses
%from. Within each cell should be an mxc matrix with number of rows, m,
%equal to the number of stimuli that you will be testing (e.g. equal to
%numel(stimuli).
%varargin - noiseCorrelationScale determining how correlated the responses
%should be across cells. The higher the number, the lower the correlation
if numel(varargin) > 0
    noiseCorrelationScale = varargin{1}; %number, the bigger the number the lower the noise correlation
    
    %simulate to calculat the noise correlation value
    rng(1);
    t1 = rand(1, runsPerStimulus);
    rng(1);
    t2 = rand(1, runsPerStimulus);
    rng('shuffle');
    t1 = t1 + (randn(1, runsPerStimulus)+0.5)*noiseCorrelationScale;
    t2 = t2 + (randn(1, runsPerStimulus)+0.5)*noiseCorrelationScale;
    
    cMat = corrcoef(t1, t2);
    noiseCorrelationCoefficient = cMat(1, 2)
    
    
else
    noiseCorrelationScale = 0;
end

trainingData = zeros(numel(stimuli)*runsPerStimulus, numel(cdfs) + 1);
for i = 1:numel(stimuli)
    trueStim = stimuli(i); %grab the true stimulus
    
    if numel(varargin) == 0
        seed = 'shuffle';
    else
        seed = randi(1000000); %reseed for each stimulus
    end
        
    for j = 1:numel(cdfs)
        c = cdfs{j};
        c = c(i, :); %grab the proper cdf from the proper cell
        %generate n random responses for each cell
        trainingData((i-1)*runsPerStimulus+1:i*runsPerStimulus, j) = randNumSpikesGenerator(c, runsPerStimulus, seed, noiseCorrelationScale);
    end
    
    trainingData((i-1)*runsPerStimulus+1:i*runsPerStimulus, numel(cdfs) + 1) = trueStim;
end
end