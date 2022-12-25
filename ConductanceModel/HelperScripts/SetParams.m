function params = SetParams(useUpCells, useDownCells)
%builds a params struct that will be used to run the capacitance model.
%User can choose whether to simulate upCells, downCells, or both. Ex: set
%useUpCells input value to 1 to use data from up cells and 0 to not use data
%from up cells. upCells and downCells cannot both be zero

if useUpCells + useDownCells == 0
    errormsg('Must include data from at least one cell type')
    return
end

%Farads - updated through 3/23/2022
upCs = 1e-12*[26.13, 77.08, 27.63, 13.69, 27.33, 52.08, 59.75, 17.8, 72.72, 67.03, 67.03, 72.04, 47.58, 12.74, 70.07, 19.37, 84.7, 69.44, 78.65, 82, 71.37, 63.7, 30.47, 57.52, 19.66, 63.7, 65.32];
downCs = 1e-12*[83.78, 28.34, 38.16, 118.6, 82.0, 89.63, 27.05, 140.1, 77.86, 111.7, 23.5, 76.32, 80.29, 85.64, 94.00, 78.65, 72.72, 74.83]; 
%upRseries = 1e6*[14.04, 16.64, 14.88, 9.945, 21.55, 37.87, 26.78, 6.843, 33.56, 25.59, 24.3, 27.95, 42.97, 56.44, 18.29, 4.681, 29.43, 28.02, 23.07, 26.42, 20.76, 31.1, 10.7, 26.58, 10.15, 29.75, 28.85];
%downRseries = 1e6*[28.31, 20.38, 11.04, 12.19, 21.96, 16.56, 7.121, 14.45, 31.8, 25.64, 7.021, 22.01, 18.3, 25.53, 24.92, 32.42, 47.76, 29.85];


[upRinput, downRinput] = GetInputResistances(); %Ohms

%Volts
upRestingVms = 1e-3*[-61.2631825000000,-53.7329043750000,-57.8621900000001,-56.4269806250000,-56.9275487500000,-46.8226906250000,-49.5415750000000,-50.7629806250000,-60.8276468750000,-53.2981343750000];
upThresholds = 1e-3*[-44.6116303080513,-42.8131431698711,-50.1307968841630,-45.1413409638600,-44.0586824775950,-36.2646492453477,-49.7336380087327,-56.0034206959113];
downRestingVms = 1e-3*[-56.5668412500001,-53.1919462499999,-51.5376512500000,-57.7276887500000,-60.3677768750000,-51.4389568750000,-51.3337675000000,-55.1412131250000,-50.9341687500000,-51.8521531250000,-60.5159781250000];
downThresholds = 1e-3*[-47.9043403256158, -44.4300303297700,-38.8380083183140,-48.2739665939588,-46.1270188066423,-45.9089223222254,-47.3296594671501,-43.5264612353869];

allCapacitances = [];
allInputResistances = [];
allRestingVms = [];
allThresholds = [];

%combine data from up cells and down cells depending on what the user
%wants.
if useUpCells
    disp('Including data from UP Cells in the model')
    allCapacitances = [allCapacitances, upCs];
    allInputResistances = [allInputResistances, upRinput];
    allRestingVms = [allRestingVms, upRestingVms];
    allThresholds = [allThresholds, upThresholds];
end
if useDownCells
    disp('Including data from DOWN Cells in the model')
    allCapacitances = [allCapacitances, downCs];
    allInputResistances = [allInputResistances, downRinput];
    allRestingVms = [allRestingVms, downRestingVms];
    allThresholds = [allThresholds, downThresholds];
end


%fill parameters for the model
params.MembraneCapacitance = median(allCapacitances);
params.InputResistance = median(allInputResistances);
params.RestingVm = median(allRestingVms);
params.SpikeThreshold = median(allThresholds);
end
