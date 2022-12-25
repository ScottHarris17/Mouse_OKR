%Processes recording epochs in data struct. Called by master_OKR_Physiology

%This script opens dataViewerApp_OKR_Physiology sequentially for each cell
%in the recording.

cells = fieldnames(data);

%open dataViewerApp for each cell
for i = size(cells, 1):-1:1
    cell_i = cells{i, 1};
    if ~strcmp(cell_i, 'centerData') %one field is centerData, don't want to look at this one here
        %% User chooses which dataviewer app to run
        h = dataViewerApp_OKR_Physiology_ObjAnalysis(data, cell_i, cellListPath); %primary analysis completed w/in dataViewerApp
        %waitfor(h) %wait for user to close app to open another instance of it for the next cell
    end
end

clear i cell_i h cells