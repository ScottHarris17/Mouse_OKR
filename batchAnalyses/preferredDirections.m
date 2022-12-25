% Builds batch analysis for OKR physiology data
load('cellList_OKR.mat') %load the batch data
tag = 'MTN'; %this is the tag that you want to identify the cells by.
cells = {}; %will hold cells that meet criteria
structs = {};

%get names of all cells that meet inclusion criteria
for i = 1:size(cellList, 1)
    s = cellList{i, 2};
    if isfield(s, 'PairedCells')
        %continue
    end
    for j = 1:numel(s.Tags)
        try
            if strcmp(s.Tags{j}, tag)
                cells{end+1} = s.cellID;
                structs{end + 1} = s;
            	break
            end
        catch
            if strcmp(s.Tags, tag)
                cells{end+1} = s.cellID;
                structs{end + 1} = s;
                break
            end
        end      
    end
end

%choose which protocols to look at
priorityOrder_1 = {'Moving_Bar'; 'Moving_Grating_Direction'}; %first field name should be one of these. Priority is in given order
priorityOrder_2 = {'_ff'; 'Extracellular'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'Extracellular', '_SP10'};
cantHave_1 = [];
cantHave_2 = {'Intracellular', 'CurrentClamp'}; %avoid anything with these tags
rigMandates = [];

noAnalyses = [];
analysesUsed = {};
pd = []; %store preferred direction
dsi = []; %store DSI
coordinates = zeros(0, 2);
cnames = {};

for i = 1:numel(structs)
    struct_i = structs{i};
    
    %grab preferred direction, and vector length
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    
    if pass
        continue
    end
    cnames{end+1} = struct_i.cellID;
    pd(end+1) = loc.Analysis_Results.PreferredDirection;
    dsi(end+1) = loc.Analysis_Results.DSI;
    try
        coordinates = [coordinates; struct_i.coordinates.polar];
    catch
        coordinates = [coordinates; [NaN NaN]];
    end
end

dsiThreshold = 0.05;
dsiAboveThreshold = dsi(dsi > dsiThreshold);

%polar plot of DSI vectors
figure()
for i = 1:numel(dsiAboveThreshold)
        polarplot([0 deg2rad(pd(i))], [0 dsiAboveThreshold(i)], 'LineWidth', 4, 'Color', [0 0 0 0.7])
        hold on
end
title('Preferred Direction of MTN labeled Cells')
hold off
%histogram of DSI magnitudes
figure()
title('Direction Selectivity Indices')
hold on
histogram(dsi, 'BinEdges', 0:0.1:1)
plot([dsiThreshold dsiThreshold], [0, 1], '--r')
hold off

%polar histogram for preferred direction
giveDegrees = 60;
availablePDs = mod(pd(dsi>dsiThreshold), 360);
upPDs = availablePDs(availablePDs > 90-giveDegrees & availablePDs < 90 + giveDegrees);
downPDs = availablePDs(availablePDs > 270 - giveDegrees & availablePDs < 270+ giveDegrees);
unwantedPDs = availablePDs(availablePDs <= 90 - giveDegrees | availablePDs >= 270 + giveDegrees |...
    (availablePDs >= 90 + giveDegrees & availablePDs <= 270 - giveDegrees));

figure()
b1 = polarhistogram(deg2rad(upPDs),'BinEdges', deg2rad([0:15:360]), 'FaceColor', 'k');
hold on
b2 = polarhistogram(deg2rad(downPDs), 'FaceColor', 'm');
b3 = polarhistogram(deg2rad(unwantedPDs), 'FaceColor', [0, 1, 1]);
b2.BinEdges = b1.BinEdges;
b3.BinEdges = b1.BinEdges;
title('Preferred Directions')
hold off

%preferred direction and location on the retina
arrowLength = 0.1;
figure()
hold on
coordinatesThreshold = coordinates(dsi>dsiThreshold, :);
for i = 1:size(coordinatesThreshold, 1)
    [x, y] = pol2cart(coordinatesThreshold(i, 1), coordinatesThreshold(i, 2));
    x = -x;
    theta = deg2rad(pd(i));
    if availablePDs(i) > 90-giveDegrees && availablePDs(i) < 90 + giveDegrees
        quiver(-x, y, arrowLength*cos(theta), arrowLength*sin(theta), '-k', 'LineWidth', 1)
    elseif availablePDs(i) > 270 - giveDegrees && availablePDs(i) < 270+ giveDegrees
        quiver(-x, y, arrowLength*cos(theta), arrowLength*sin(theta), '-m', 'LineWidth', 1)
    else
        quiver(-x, y, arrowLength*cos(theta), arrowLength*sin(theta), 'Color', [0, 1, 1], 'LineWidth', 1)
    end
end
circleDraw(0, 0, 1);
title('Preferred Direction and Soma Location')
hold off

function h = circleDraw(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, '-k');
scatter(x, y, '+')
hold off
end
                

        