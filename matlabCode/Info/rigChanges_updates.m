%make a structure called rig. Update the rig subfields (attributes)
%according to changes that you make to the rig. Each subfield is a
%particular trait of the rig, for example, the diameter of the LED. Place
%as second order subfields the date that the trait was changed and the
%value that it was changed to. Update and rerun this script every time you
%make a change to the rig.
rig = struct();

rig.experimentStage.d20201215 = 1; %start of recordings from MTN labeled cells
rig.experimentStage.d20210725 = 2; %End of cell attached recordings and start of whole cell recordings
rig.experimentStage.d20211201 = 3; %Problems with cell attached recordings between 11/1/2021 and 12/1/2021
rig.LEDSpotDiameter.d20200901 = 500; %um
rig.LEDSpotDiameter.d20210817 = 300; %um






save('rigChanges.mat', 'rig'); %save the structure as a .mat file. This file can then be used in batch analysis to seperate out data that was run on different rig configurations if desired