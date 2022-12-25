% add pipette solution to those cells that don't have it in their meta
% data.
for i = 1:size(downs, 1)
    cname = downs{i, 1};
    for j = 1:size(cellList, 1)
        name_j = cellList{j, 1};
        if strcmp(cname, name_j)
            break
        end
    end
    cellListIndex = j;
    
    struct_temp = downs{i, 2};
    fLocation = struct_temp.fileLocation;
    fname = [fLocation , cname, '_Analysis.mat'];
    load(fname)
    tags = CellParameters.Tags;
    
    found = 0;
    for j = 1:numel(tags)
        tag_j = tags{j};
        if strcmp('DOWN Cell', tag_j)
            found = 1;
            break
        end
    end
    if found
        continue
    end
    
    tags = [tags, {'DOWN Cell'}];
    CellParameters.Tags = tags;
    
    cellList{cellListIndex, 2} = CellParameters;
    save(fname, 'CellParameters')
    clear CellParameters fLocation struct_temp tags j cellListIndex cname
end
