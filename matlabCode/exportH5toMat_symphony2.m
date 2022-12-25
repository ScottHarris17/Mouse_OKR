%This script automates the job of exporting .h5 files to a clarinet export.
%exports are saved as .mat files with the cell name + extension
%'_ClarinetExport' to the same folder as the .h5 file. This script uses code
%that is installed with clarinet.

%20221202 - SH
%% User selects .h5 file
function exportH5toMat_symphony2(fName, PathName)
cd(PathName);
fName = [fName, '.h5'];
%% import h5 file
cells = [];
import parsers.*
version = SymphonyParser.getVersion(fName);
if version == 2
    ref = SymphonyV2Parser(fName);
else
    ref = SymphonyV1Parser(fName);
end
ref.parse;
data = ref.getResult;

cells = [];
for j = 1:numel(data)
  if isempty(cells)
      cells = data{j};
  else
      cells(end+1) = data{j};
  end
end

%% save a .mat file for each cell in the h5 file
for j = 1:numel(cells)
    disp(['Saving data for cell # ', num2str(j)])
    
    FileName = cells(j).get('h5File');
    CellName = cells(j).get('label');
    if iscell(CellName)
        CellName = [CellName{:}];
    end
    for i = 1:numel(cells(j).epochs)
        epoch = cells(j).epochs(i);
        meta = epoch.toStructure;
        device = epoch.get('devices'); %change this if you ever use anything other than Amp_Ch1
        response = epoch.getResponse(device{1});
        epochs(i).epoch = response.quantity'; % Retrieve epoch data
        epochs(i).meta  = meta;
    end
    retinaName = FileName(strfind(FileName, '_')+1:strfind(FileName, '.')-1);
    cname = [retinaName 'c' CellName(2:end)];
    savename = [cname '_ClarinetExport'];
    save(savename, 'epochs');
    clear epochs
end                 
end


