%example of changing the pipette solution for a number of cells based off
%of their completed analysis names
for i = 210:size(cellList, 1)
    cname = cellList{i, 1};
    struct_temp = cellList{i, 2};
    fLocation = struct_temp.fileLocation;
    fname = [fLocation , cname, '_Analysis.mat'];
    load(fname)
    analyses = fieldnames(CellParameters.Analyses_Completed);
    for j = 1:numel(analyses)
        analysis_j = analyses{j};
        subFields = fieldnames(CellParameters.Analyses_Completed.(analysis_j));
        for k = 1:numel(subFields)
            subfield_k = subFields{k};
            if contains(subfield_k, 'Cesium')
                disp([cname, ' ', analysis_j, ' ', subfield_k ])
                CellParameters.Analyses_Completed.(analysis_j).(subfield_k).meta.pipetteSolution = 'Cesium';
                CellParameters.Analyses_Completed.(analysis_j).(subfield_k).Analysis_Results.meta.pipetteSolution = 'Cesium';
            end
        end
    end
    %save(fname, 'CellParameters')
    cellList{i, 2} = CellParameters;
end
    
%     cellList{cellListIndex, 2} = CellParameters;
%     save(fname, 'CellParameters')
% end
