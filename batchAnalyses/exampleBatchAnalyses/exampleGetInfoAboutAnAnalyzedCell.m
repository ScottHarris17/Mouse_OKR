load('C:\Users\mrsco\Box\DataAndAnalysis\labData\OKR\batch\cellList_OKR.mat')
cellID = 'SHOKR76Lc12'
excitation = 0;

load('cellList_OKR.mat')
for i = 1:size(cellList, 1)
    cname = cellList{i, 1};
    if strcmp(cname, cellID)
        struct = cellList{i, 2};
        break
    end
end


%choose which protocols to look at
priorityOrder_1 = {'Moving_Bar'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; 'Intracellular'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
if excitation
    mustHave_2 = {'Intracellular', 'Excitation'};
else
    mustHave_2 = {'Intracellular', 'Inhibition'};
end
cantHave_1 = [];
cantHave_2 = {'Extracellular'}; %avoid anything with these tags
rigMandates = []; %use only cells with LED spot diameter = 500

[loc, pass] = getLoc(struct, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);

disp(['pass = ', num2str(pass)])

if ~pass
    names = fieldnames(struct.Analyses_Completed.Moving_Bar);
    disp('names = '); disp(names)
    
    epochs = loc.Analysis_Results.EpochNumbers;
    avg = loc.Analysis_Results.meanTracesByOrientation.meanTraces';
    
    figure()
    hold on
    title(cellID)
    plot(avg);
    legend({'0', '315', '270', '225', '180', '135', '90', '45'})
    
    disp(['epochNums = ', num2str(epochs)]);
end
