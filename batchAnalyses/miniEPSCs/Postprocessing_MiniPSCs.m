% Adopted from PCA analysis written by Luca Della Santina
% Scott Harris 10/23/2021

% Postprocesses minis taken from igor export text file with file name
% _cutouts.txt
%
% input: Matrix of waveforms (rows=different traces cols=timepoints)
%        ** IMPORTANT ** Baseline of the imported waveforms must be zero 
%                        otherwise amplitudes are not calculated correctly.
%
% output : wavesPCA
%           |- wavesSelected: contains user-selected waves
%           |- wavesRejected: contains all rejected waves
%           |- settings: settings used for analysis and plotting
%           |- freqSelected: frequency of selected waves in the recording
%           |- ampSelectedAvg: mean amplitude of selected waves
%           |- ampSelectedSD: St.dev. of selected waves' amplitudes
%
% -------------------------------------------------------------------------
% Version 4.1                               2017-06-26  Luca Della Santina
%
% + New plotting modes: PC1 vs PC2, and PC1 vs Time-to-peak (negative pk)
% + Allows to select wavesPCA results as input to refine previous analysis
% - Removed default k-means clustering, user can select it from menu
%
% Version 4.0                               2017-06-23  Luca Della Santina 
%
% + New selection mode: "invalidate data", allows recursive invalidation
% + New button redo k-means clustering on remaining valid data
%
% Version 3.0                               2017-06-23  Luca Della Santina 
%
% + When mouse cursor move, display closest waveform in right-top panel
% + When clicking on point, store closest waveform in right-bottom panel
% + Added gui buttons for manually selecting waveforms and exit
% % Fixed: rejected waveform were plotted with SEM instead of SD
%
% Version 2.0                               2017-06-21  Luca Della Santina 
%
% + Plots average waveform for selected and rejected waves with SEM
% + Plots average amplitude and frequency (control vs treat experiments)
% + Asks user to insert time scale and total recording time
% + Calculate and save mean amplitude of selected waveform (negative peak)
% + Calculate frequency of selected waveforms
% % Allows user to refine k-mean cluster by drawing polygon
%
% Version 1.0                               2017-06-19  Luca Della Santina 
%
% + Plots k-means clustering next to manual clustering
% -------------------------------------------------------------------------

% Use this line to load Igor-exported waves .txt

%input cell names you want to analyze
analyzeTheseCells = {'SHOKR85Lc11',};


cellListPath = 'C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\batch\cellList_OKR.mat';
load(cellListPath)

%%
for i = 1:numel(analyzeTheseCells)
    cell_i = analyzeTheseCells{i};
    disp(['Beginning Post processing for cell ' cell_i])

    location = strcmp(cellList(:, 1), cell_i);
    struct_i = cellList{location, 2};
    
    fileLocation = struct_i.fileLocation;
    igorExportPath = [fileLocation, '\MiniPSCs'];
    try
        docs = dir(igorExportPath);
    catch
        warning(['Cannot find MiniPSCs directory for cell' cell_i])
    end
    
    availableIgorExports = {};
    analysisNames = {};
    for j = 1:numel(docs)
        name_j = docs(j).name;
        if contains(name_j, '_cutout.txt')
            bookendIndex = strfind(name_j, '_cutout.txt');
            analysisNames{1, end + 1} = name_j(1:bookendIndex-1);
            availableIgorExports{1, end + 1} = fullfile(igorExportPath, name_j);            
        end
    end
    
    recordingLength = 1; %s
    sampleRate = 10; %sampling rat (kHz);
    
    %cycle through each available Igor Export for this cell and do the
    %manual PCA clustering.
    for j = 1:numel(availableIgorExports)
        %add waves
        currentAnalysis = analysisNames{j};
        disp(['Cell ' cell_i ': ' currentAnalysis])

        fname = availableIgorExports{j};
        waves=dlmread(fname)';
        tmpWaves  = waves;
        tmpSettings.totalLen = recordingLength;
        tmpSettings.samplingRate = sampleRate;
        [coef, score, latent] = pca(tmpWaves); % newer function
        score = score(:,1:2); % just keep the score for the first 2 PCA components
        tmpOriginalScore = score;
        tmpCumSum = cumsum(latent)./sum(latent) * 100;
        tmpH = figure('Name', ['PCA analysis - First two PCA components accounting for '  num2str(tmpCumSum(2),2) '% of variance']);
        set(tmpH, 'Position', [100 200 1200 500]);
        
        
        % Gui:
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'PC1 vs PC2','Units','normalized',...
            'Position',[0.01 0.85 0.08 0.08],'Visible','on',...
            'CallBack',['redoPC1vsPC2=1; uiresume']);
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'PC1 vs Peak time','Units','normalized',...
            'Position',[0.01 0.75 0.08 0.08],'Visible','on',...
            'CallBack',['redoPC1vsPeakTime=1; uiresume']);
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'Select Data','Units','normalized',...
            'Position',[0.01 0.65 0.08 0.08],'Visible','on',...
            'CallBack',['selManual=1; uiresume']);
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'Invalidate Data','Units','normalized',...
            'Position',[0.01 0.55 0.08 0.08],'Visible','on',...
            'CallBack',['invManual=1; uiresume']);
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'k-means','Units','normalized',...
            'Position',[0.01 0.45 0.08 0.08],'Visible','on',...
            'CallBack',['redoKmeans=1; uiresume']);
        uicontrol('Parent',tmpH,'Style','pushbutton','String', 'Exit','Units','normalized',...
            'Position',[0.01 0.09 0.08 0.08],'Visible','on',...
            'CallBack',['exitFlag=1; uiresume']);
        set (gcf, 'WindowButtonMotionFcn', ['mouseMove=1; uiresume']);
        set (gcf, 'WindowButtonDownFcn', ['mouseDown=1; uiresume']);
        
        tmpScore = score; % tmpScore is a temporary list of valid waveforms, allows refining by manually invalidating
        exitFlag=0;
        selManual=0;
        invManual=0;
        redoKmeans=0;
        redoPC1vsPC2=1;
        redoPC1vsPeakTime=0;
        mouseMove=0;
        mouseDown=0;
        
        while exitFlag >= 0
            if exitFlag
                % Close the program and clear the workspace
                clear tmp* ans selManual invManual redo* exitFlag mouse*;
                clear latent pcvars score X Y zscore in out hpoly C idx coef;
                close all;
                break
            elseif redoPC1vsPC2
                % Recluster current valid data with k-means algorithm
                score = tmpOriginalScore; % reload score from original PCA results
                tmpScore = score;
                subplot(2,2,[1 3]);
                hold off; scatter(score(:,1),score(:,2),140,'k.'); hold on;
                xlabel(['PCA component 1 - (' num2str(round(latent(1)/sum(latent)*100)) '% of variance)']);
                ylabel(['PCA component 2 - (' num2str(round(latent(2)/sum(latent)*100)) '% of variance)']);
                set(gca,'TickDir','out');
                tmpGCA = gca;
                redoPC1vsPC2=0;
            elseif redoPC1vsPeakTime
                % Recluster current valid data with k-means algorithm
                score = tmpOriginalScore; % reload score from original PCA results
                [tmpPeaksAmp tmpPeaksTime] =  min(tmpWaves, [], 2);
                tmpPeaksTime = tmpPeaksTime / tmpSettings.samplingRate;
                score = [score(:,1) tmpPeaksTime];
                tmpScore = score;
                subplot(2,2,[1 3]);
                hold off; scatter(score(:,1),score(:,2),140,'k.'); hold on;
                xlabel(['PCA component 1 - (' num2str(round(latent(1)/sum(latent)*100)) '% of variance)']);
                ylabel('Time to negative peak');
                set(gca,'TickDir','out');
                tmpGCA = gca;
                redoPC1vsPeakTime=0;
            elseif redoKmeans
                % Recluster current valid data with k-means algorithm
                subplot(2,2,[1 3]);
                [idx,C] = kmeans(tmpWaves(find(ismember(score(:,1), tmpScore(:,1))),:), 2);
                scatter(tmpScore(find(idx==1),1),tmpScore(find(idx==1),2),140,'r.'); hold on;
                scatter(tmpScore(find(idx==2),1),tmpScore(find(idx==2),2),140,'b.'); hold on;
                % wavesKmeans1 = tmpWaves(find(idx==1),:); % Save k-means results
                % wavesKmeans2 = tmpWaves(find(idx==2),:);
                redoKmeans=0;
            elseif invManual
                % Refine automatic cluster by drawing polygon around data to invalidate
                subplot(2,2,[1 3]);
                title('Refine k-means clustering by drawing polygon around invalid data');
                [X,Y]=getline(gcf,'closed'); % User selects data by drawing polygon
                hpoly=plot(X,Y,'b');
                in=find(inpolygon(tmpScore(:,1),tmpScore(:,2),X,Y)); % Find points within polygon
                out=setdiff(1:size(tmpScore,1), in)'; % Find points outside of polygon
                tmpScore = tmpScore(out,:);
                hold off; scatter(tmpScore(:,1),tmpScore(:,2),140,'k.'); hold on;
                pause(1);
                set (gcf, 'WindowButtonMotionFcn', ['mouseMove=1; uiresume']);
                set (gcf, 'WindowButtonDownFcn', ['mouseDown=1; uiresume']);
                invManual=0;
            elseif selManual
                % Manually refine automatic cluster by drawing polygon around valid data
                subplot(2,2,[1 3]);
                title('Refine k-means clustering by drawing polygon around valid data');
                [X,Y]=getline(gcf,'closed'); % User selects data by drawing polygon
                
                hpoly=plot(X,Y,'b');
                in=find(inpolygon(tmpScore(:,1),tmpScore(:,2),X,Y)); % Find points within polygon within the refined tmpScore
                out=setdiff(1:size(score,1), in)'; % Find points outside of polygon from the entire list contained in score
                scatter(tmpScore(in,1),tmpScore(in,2),140,'go');
                legend('k-means group #1', 'k-means group #2', 'User selected',...
                    'Location','best');
                %scatter(score(out,1),score(out,2),140,'r.');
                pause(1);
                
                % Plot average waveform of selected and rejected waves
                tmpSelected = tmpWaves(in,:);
                tmpRejected = tmpWaves(out,:);
                
                subplot(2,2,2);
                tmpMean=mean(tmpSelected,1);
                tmpSD=std(tmpSelected,1);
                tmpSEM=tmpSD/sqrt(size(tmpSelected,1));
                tmpTimeX = 1:size(tmpWaves,2); % Scale time axis for sampling rate
                tmpTimeX = tmpTimeX / tmpSettings.samplingRate;
                plot(tmpTimeX, tmpMean, 'b', 'LineWidth',1); hold on;
                plot(tmpTimeX, tmpMean+tmpSD, '--', 'Color',[0.5,0.5,0.5]);
                plot(tmpTimeX, tmpMean-tmpSD, '--', 'Color',[0.5,0.5,0.5]);
                set(gca,'box','off');
                xlim([0 max(tmpTimeX)]);
                ylabel('Current (pA)'); title('Selected waveforms');
                legend('Average','Standard Deviation','Location','southeast');
                
                subplot(2,2,4);
                tmpMean=mean(tmpRejected,1);
                tmpSD=std(tmpRejected,1);
                tmpSEM=tmpSD/sqrt(size(tmpRejected,1));
                plot(tmpTimeX, tmpMean, 'b', 'LineWidth',1); hold on;
                plot(tmpTimeX, tmpMean+tmpSD, '--', 'Color',[0.5,0.5,0.5]);
                plot(tmpTimeX, tmpMean-tmpSD, '--', 'Color',[0.5,0.5,0.5]);
                set(gca,'box','off');
                xlim([0 max(tmpTimeX)]);
                xlabel('Time (ms)'); ylabel('Current (pA)'); title('Rejected waveforms');
                legend('Average','Standard Deviation','Location','southeast');
                
                % Save analysis results into wavesPCA structure
                wavesPCA.wavesSelected  = tmpSelected;
                wavesPCA.wavesRejected  = tmpRejected;
                wavesPCA.settings       = tmpSettings;
                wavesPCA.freqSelected   = size(tmpSelected,1)/tmpSettings.totalLen;
                wavesPCA.ampSelectedAvg = abs(mean(min(tmpSelected, [], 2)));
                wavesPCA.ampSelectedSD  = std(min(tmpSelected, [], 2));
                
                % release mouse events detection, otherwise strange gui behavior
                set (gcf, 'WindowButtonMotionFcn', '');
                set (gcf, 'WindowButtonDownFcn', '');
            elseif mouseMove
                % Detect mouse position on the plot (using plot coordinates)
                mousePos = get (tmpGCA, 'CurrentPoint'); mousePos = mousePos(1,1:2);
                % Compute Euclidean distances of scores to the mouse pointer:
                tmpDist = sqrt(sum(bsxfun(@minus, tmpScore, mousePos).^2,2));
                % Find the smallest distance and use that as an index into B:
                tmpClosestScore = tmpScore(tmpDist==min(tmpDist),:);
                % Plot the closest waveform
                subplot(2,2,2);
                tmpClosestWave = tmpWaves(find(tmpScore(:,1) == tmpClosestScore(:,1)),:);
                tmpTimeX = 1:size(tmpWaves,2); % Scale time axis for sampling rate
                tmpTimeX = tmpTimeX / tmpSettings.samplingRate;
                plot(tmpTimeX, tmpClosestWave, 'b', 'LineWidth',1);
                xlim([0 max(tmpTimeX)]);
                title(['Closest waveform [X, Y] to mouse position (X,Y) = [',...
                    num2str(tmpClosestScore(1,1),3), ', ',num2str(tmpClosestScore(1,2),3), '] (',...
                    num2str(mousePos(1,1),3), ', ', num2str(mousePos(1,2),3), ')']);
                ylabel('Current (pA)');
                xlabel('Time (ms)');
                mouseMove=0;
            elseif mouseDown
                % Detect mouse position on the plot (using plot coordinates)
                mousePos = get (tmpGCA, 'CurrentPoint'); mousePos = mousePos(1,1:2);
                % Compute Euclidean distances of scores to the mouse pointer:
                tmpDist = sqrt(sum(bsxfun(@minus, tmpScore, mousePos).^2,2));
                % Find the smallest distance and use that as an index into B:
                tmpClosestScore = tmpScore(tmpDist==min(tmpDist),:);
                % Plot the closest waveform
                subplot(2,2,4);
                tmpClosestWave = tmpWaves(find(tmpScore(:,1) == tmpClosestScore(:,1)),:);
                tmpTimeX = 1:size(tmpWaves,2); % Scale time axis for sampling rate
                tmpTimeX = tmpTimeX / tmpSettings.samplingRate;
                plot(tmpTimeX, tmpClosestWave, 'r', 'LineWidth',1);
                xlim([0 max(tmpTimeX)]);
                title(['Last clicked waveform [X, Y] close to mouse position (X,Y) = [',...
                    num2str(tmpClosestScore(1,1),3), ', ',num2str(tmpClosestScore(1,2),3), '] (',...
                    num2str(mousePos(1,1),3), ', ', num2str(mousePos(1,2),3), ')']);
                ylabel('Current (pA)');
                xlabel('Time (ms)');
                mouseDown=0;
            end
            uiwait;
        end
        

        %check that the field exists
        
        if ~isfield(struct_i.Analyses_Completed.MiniPSCs, (currentAnalysis))
            error(['Something is wrong. The ' currentAnalysis ' for cell ' cell_i ' could not be located in the structure from cellList'])
        end
        loc = struct_i.Analyses_Completed.MiniPSCs.(currentAnalysis);
        
        %grab the peak value for every detected event
        normedWaves = wavesPCA.wavesSelected;
        for i = 1:size(normedWaves, 1)
            normedWaves(i, :) = normedWaves(i, :)-mean(normedWaves(i, 1:10));
        end
        
        allPeaks = max(abs(normedWaves), [], 2); %always normalized to be positive
        
        loc.Analysis_Results.numWaves = wavesPCA.freqSelected;
        loc.Analysis_Results.numWavesRejected = size(wavesPCA.wavesRejected, 1);
        loc.Analysis_Results.allPeaks = allPeaks; %
        loc.Analysis_Results.meanWave = mean(wavesPCA.wavesSelected);
        loc.Analysis_Results.meanWaveAmplitude = wavesPCA.ampSelectedAvg;
        loc.Analysis_Results.stdWaveAmplitude = wavesPCA.ampSelectedSD;
        loc.Analysis_Results.postprocesseingComplete = 1;
        struct_i.Analyses_Completed.MiniPSCs.(currentAnalysis) = loc;
    end
    %once done with all the analyses for this cell, put everything back
    %where it belongs and save
    cellList{location, 2} = struct_i;
    
    %save the updated analysis structure back to its proper location
    CellParameters = struct_i;
    save([struct_i.fileLocation, cell_i, '_Analysis.mat'], 'CellParameters')
    disp(['Saved all data for ' cell_i])
end

%% Save the cell list
save(cellListPath, 'cellList')
disp('Saved cellList.mat')
