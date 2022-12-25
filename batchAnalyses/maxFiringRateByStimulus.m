%% Analyze the firing rates of each cell during grating and flash stimuli

% This analysis  It compares the average maximum firing rate during the
% flash to the average maximum firing rate during the moving grating 
%(in the preferred direction).

cellListFName = "C:/Users/mrsco/Box/DataAndAnalysis/labData/OKR/batch/cellList_OKR.mat";
 
windowWidth = 50; %in ms, how wide will the sliding window be overwhich spikes are counted together?
stepSize = 10; %in ms, how big will the step size be that the window moves by on each iteration?
DSI_cutoff = 0.1; %minimum DSI to count as direction selective

%%
load(cellListFName); %loads cellList to workspace

%initialize arrays to hold results for each cell. These will be 2xN
%matrices where the first column is the average maximum firing rate for
%each cell in response to the flash and the second column is the average
%maximum firing rate for each cell in response to the moving grating.
ONcells_firingRate = [];
ONcells_isi = {};
ON_DS = [];
OFFcells_firingRate = [];
OFFcells_isi = {};
OFF_DS = [];
ON_OFFcells_firingRate = [];
ON_OFFcells_isi = {};
ON_OFF_DS = [];


for i = 1:size(cellList, 1) %iterate for each cell
    cp = cellList{i, 2}; %cp is the cellparameters structure for the current cell
    flashAnalyses = cp.Analyses_Completed.Flash;
    
    try
        movingGratingAnalyses = cp.Analyses_Completed.Moving_Grating_Direction;
    catch
        warning(['Not analyzing data for ', cp.cellID , ' - no Moving_Grating_Direction analysis completed for this cell']);
    end
    
    flashTypes = fieldnames(flashAnalyses); %returns an nx1 cell array of all the analyses that have been done for flash
    movingGratingTypes = fieldnames(movingGratingAnalyses);
    
    %figure out which flash analyses to use. Priority is given to a flash
    %analysis with the suffix '_ff', and then to a flash analysis with the
    %suffix '_all'.
    flashFieldName = '';
    for j = 1:numel(flashTypes)
        if contains(flashTypes{j}, '_ff')
            flashFieldName = flashTypes{j};
        elseif strfind(flashTypes{j}, '_all') && strcmp(flashFieldName, '')
            flashFieldName = flashTypes{j};
        end
    end
    
    %figure out which moving grating analysis to use. Priority is given to
    %a stimulus with the suffix '_ff'. If this doesn't exist, the cell
    %isn't analyzed.
    gratingFieldName = '';
    for j = 1:numel(movingGratingTypes)
        if contains(movingGratingTypes{j}, '_ff')
            gratingFieldName = movingGratingTypes{j};
        end
    end
    
    if strcmp(gratingFieldName, '')
        warning(['Analysis not performed for ', cellList{i, 1}, ' No final moving grating direction stimulus found'])
        continue %jump to next cell in cellList
    end
    
    %grab just the analysis results you want
    flashAnalysis = flashAnalyses.(flashFieldName);
    gratingAnalysis = movingGratingAnalyses.(gratingFieldName);
    
    %grab relevant data for flash and grating stimuli (it's possible that
    %not all flashes will have the same background and stim light levels.
    %For now, analyze them all together. In the future, maybe be careful in
    %the dataViewerApp about grouping certain ones together based on these
    %parameters).
    preTimeFlash = flashAnalysis.meta.preTime;
    stimTimeFlash = flashAnalysis.meta.stimTime;
    tailTimeFlash = flashAnalysis.meta.tailTime;
    spikeTimesFlash = flashAnalysis.Analysis_Results.allSpikeTimes;
    
    preTimeGrating = gratingAnalysis.meta.preTime;
    stimTimeGrating = gratingAnalysis.meta.stimTime;
    tailTimeGrating = gratingAnalysis.meta.tailTime;
    spikeTimesGrating = gratingAnalysis.Analysis_Results.allSpikeTimes;
   
    %grab average maximum spike rate for this cell to the moving grating stimulus
    %in the preferred direction.
    preferredDirection = gratingAnalysis.Analysis_Results.PreferredDirection;
    if preferredDirection < 0
        preferredDirection = 360 + preferredDirection;
    end
    
    %find the orienation used in the experiment that is closest to the
    %preferred direction for the cell.
%     orientations = gratingAnalysis.meta.orientation; %row vector of all orientations used
%     orientationDifs = orientations - preferredDirection;
%     [~, I_min] = min(abs(orientationDifs));
%     bestOrientation = orientations(I_min); 
    
    [~, ind] = max(gratingAnalysis.Analysis_Results.mean_spikesByOrientation);
    bestOrientation = gratingAnalysis.meta.orientation(ind);
    
    epochsWithBestOrientation = find(gratingAnalysis.Analysis_Results.orientationByEpoch == bestOrientation);
    
    %grab the spikes for the desired epochs
    gratingMaxRates = [];
    grating_isi = [];
    for j = 1:numel(epochsWithBestOrientation)
        spikes = gratingAnalysis.Analysis_Results.allSpikeTimes{epochsWithBestOrientation(j)};
        relevantSpikes = (spikes > (preTimeGrating + 50)) &...
                (spikes < (preTimeGrating+stimTimeGrating+ 50)); 
        
        relevantSpikeTimes = spikes(relevantSpikes);
        gratingMaxRates = [gratingMaxRates, maxSpikeRate(relevantSpikeTimes, windowWidth, stepSize)];
        grating_isi = [grating_isi, diff(relevantSpikeTimes)];
    end
    avgGratingMaxRates = mean(gratingMaxRates);
    
    DSI_classification = gratingAnalysis.Analysis_Results.DSI > DSI_cutoff; %takes 1 if cell is DS, 0 if cell is not
    
    %get max spike rates for flashes as well, then record average maxes for
    %flash and grating intot the relevant matrices initialized just before
    %the outter-most for-loop.
    if strcmp(cp.Cell_Type.Physiology, 'ON')
        allMaxRatesFlash = [];
        flash_isi = [];
        for j = 1:numel(spikeTimesFlash)
            spikes = spikeTimesFlash{j}; %time of each spike in ms
            relevantSpikes = (spikes > (preTimeFlash + 50)) &...
                (spikes < (preTimeFlash+stimTimeFlash+ 50));
            
            relevantSpikeTimes = spikes(relevantSpikes);% grab only the relevant spikes
            allMaxRatesFlash = [allMaxRatesFlash maxSpikeRate(relevantSpikeTimes, windowWidth, stepSize)];
            flash_isi = [flash_isi, diff(relevantSpikeTimes)];
        end
        avgMaxRateFlash = mean(allMaxRatesFlash);
        ONcells_firingRate = [ONcells_firingRate;avgMaxRateFlash, avgGratingMaxRates];
        ONcells_isi(end+1, 1) = {flash_isi};
        ONcells_isi(end, 2) = {grating_isi};
        ON_DS = [ON_DS, DSI_classification];
        
    elseif strcmp(cp.Cell_Type.Physiology, 'OFF')
        allMaxRatesFlash = [];
        flash_isi = [];
        for j = 1:numel(spikeTimesFlash)
            spikes = spikeTimesFlash{j}; %time of each spike in ms
            relevantSpikes = (spikes > (preTimeFlash + stimTimeFlash + 25));
            
            relevantSpikeTimes = spikes(relevantSpikes);% grab only the relevant spikes
            allMaxRatesFlash = [allMaxRatesFlash maxSpikeRate(relevantSpikeTimes, windowWidth, stepSize)];
            flash_isi = [flash_isi, diff(relevantSpikeTimes)];
        end
        avgMaxRateFlash = mean(allMaxRatesFlash);
        OFFcells_firingRate = [OFFcells_firingRate;avgMaxRateFlash, avgGratingMaxRates];
        OFFcells_isi(end+1, 1) = {flash_isi};
        OFFcells_isi(end, 2) = {grating_isi};
        OFF_DS = [OFF_DS, DSI_classification];
        
        
    elseif strcmp(cp.Cell_Type.Physiology, 'ON-OFF')
        allMaxRatesFlash = [];
        flash_isi = [];
        for j = 1:numel(spikeTimesFlash)
            spikes = spikeTimesFlash{j}; %time of each spike in ms
            relevantSpikes = (spikes > (preTimeFlash + 25));
            
            relevantSpikeTimes = spikes(relevantSpikes);% grab only the relevant spikes
            allMaxRatesFlash = [allMaxRatesFlash maxSpikeRate(relevantSpikeTimes, windowWidth, stepSize)];
            flash_isi = [flash_isi, diff(relevantSpikeTimes)];
        end
        avgMaxRateFlash = mean(allMaxRatesFlash);
        ON_OFFcells_firingRate = [ON_OFFcells_firingRate;avgMaxRateFlash, avgGratingMaxRates];
        ON_OFFcells_isi(end+1, 1) = {flash_isi};
        ON_OFFcells_isi(end, 2) = {grating_isi};
        ON_OFF_DS = [ON_OFF_DS, DSI_classification];
        
        
        
    else
        warning(['No cell type listed for ' ,cellList{i, 1}, ' Analysis not performed'])
        continue
    end
end
clearvars -except cellList OFFcells_firingRate OFFcells_isi...
    ON_OFFcells_firingRate ON_OFFcells_isi ONcells_firingRate ONcells_isi...
    ON_DS OFF_DS ON_OFF_DS

%% build max spike rate figures for ON, ON-OFF, and OFF cells
for i = 1:3 %loop once for each cell type
    figure
    if i == 1
        title('Average Maximum Firing Rate - ON Cells')
        hold on
        toPlot = ONcells_firingRate;
        DS_cells = ON_DS;
        pON = signrank(ONcells_firingRate(:, 1), ONcells_firingRate(:, 2))
    elseif i == 2
        title('Average Maximum Firing Rate - OFF cells')
        hold on
        toPlot = OFFcells_firingRate;
        DS_cells = OFF_DS;
        pOFF = signrank(OFFcells_firingRate(:, 1), OFFcells_firingRate(:, 2))
    else
        title('Average Maximum Firing Rate - ON-OFF cells')
        hold on
        toPlot = ON_OFFcells_firingRate;
        DS_cells = ON_OFF_DS;
        pON_OFF = signrank(ON_OFFcells_firingRate(:, 1), ON_OFFcells_firingRate(:, 2))
    end
    
    for j = 1:size(toPlot, 1)
        if DS_cells(j)
            color = 'r';
        else
            color = 'k';
        end
        p = plot([1, 2], [toPlot(j, 1), toPlot(j, 2)],...
             color, 'Marker','.', 'MarkerFaceColor', color,...
             'MarkerEdgeColor', color, 'MarkerSize', 20, 'LineWidth', 2);
    end
    xlim([0.5 2.5])
    xticks([1, 2])
    xticklabels({'Flash', 'Grating'})
    ylabel('Hz')
    hold off
end
    
%% build interspike interval figures
% for i = 1:3
%     if i == 1
%         toPlot = ONcells_isi;
%         titleString = 'ON cell';
%     elseif i == 2
%         toPlot = OFFcells_isi;
%         titleString = 'OFF cell';
%     elseif i == 3
%         toPlot = ON_OFFcells_isi;
%         titleString = 'ON-OFF cell';
%     end
%     
%     for j = 1:size(toPlot, 1)
%         flash_isi_distribution = toPlot{j, 1};
%         grating_isi_distribution = toPlot{j, 2};
%         
%         figure
%         histogram(flash_isi_distribution,'FaceColor', 'k');
%         hold on
%         histogram(grating_isi_distribution, 'FaceColor', 'r', 'FaceAlpha', 0.6);
%         title(titleString)
%         hold off
%     end
% end
