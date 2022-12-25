%% set gains and holding potentials
orientations =   [0 45 90 135 180 225 270 315];
ExcitationGain = [1  1  1  1   1   1   1   1 ]; %gain for each stimulus direction
ExcitationGainWeight = 1; %multiplied by the gains in all directions
InhibitionGain = [1  1  1  1   1   1   1   1 ];
InhibitionGainWeight = 1;
%Choose to use data from UP Cells, DOWN Cells or both (1 = use data, 0 =
%don't use data).
useUpCells = 1;
useDownCells = 1;

ExcRev = 7; %mV
InhRev = -63; %mV
ExcHold = -63; %mV
InhHold = 7; %mV %recordings were made at +10 as indicated on multiclamp

timestep = 1e-3; %time in seconds between each data point in the conductance traces you'll use (this is the integration time step)

%% Set initial parameters based on emperical data
params = SetParams(useUpCells, useDownCells);
params.ExcitationGain = ExcitationGain;
params.ExcitationGainWeight = ExcitationGainWeight;
params.InhibitionGain = InhibitionGain;
params.InhibitionGainWeight = InhibitionGainWeight;
params.ExcitationReversal = ExcRev*1e-3; %convert to Volts
params.InhibitionReversal = InhRev*1e-3; %convert to Volts
params.ExcitationHold = ExcHold*1e-3; %convert to Volts
params.InhibitionHold = InhHold*1e-3; %convert to Volts
params.orientations = orientations;

%% grab traces for excitation and inhibition
[excitationMeans, inhibitionMeans] = GetVClampTraces(useUpCells, useDownCells);
%turn into conductances and change from pA to Amps
params.ExcitationConductances = 1e-12*excitationMeans./(params.ExcitationHold - params.ExcitationReversal);
params.InhibitionConductances = 1e-12*inhibitionMeans./(params.InhibitionHold - params.InhibitionReversal);

%% Use params to generate membrane dv/dt at each time point
ExcitationGainWeights = [0.25:0.01:3];
DSIs = zeros(2, numel(ExcitationGainWeights));
Areas = zeros(2, numel(ExcitationGainWeights));
NormalizedAreas = zeros(2, numel(ExcitationGainWeights));
VonMisesMeans = zeros(2, numel(ExcitationGainWeights));
VonMisesWidths = zeros(2, numel(ExcitationGainWeights));
VonMisesErrors = zeros(2, numel(ExcitationGainWeights));
PDs = zeros(2, numel(ExcitationGainWeights));
PDResponses = zeros(2, numel(ExcitationGainWeights));
NullResponses = zeros(2, numel(ExcitationGainWeights));

params.SpikeThreshold = -0.0462;
%plot all conductances at gain 1
% figure
% hold on
% title('Conductances')
% for i = 1:8
%     plot(params.InhibitionConductances(i, :))
% end
% plot(params.ExcitationConductances(1, :), 'g')
% legend({'0', '45', '90', '135', '180', '225', '270', '315', 'Excitation'})

g = waitbar(0,'Simulation in progress');

for i = 1:numel(ExcitationGainWeights)
    waitbar(i/numel(ExcitationGainWeights), g, 'Simulation in progress')
    params.ExcitationGainWeight = ExcitationGainWeights(i);
    %params.InhibitionGainWeight = ExcitationGainWeights(i);
    outputs = RunConductanceSimulation(params, timestep);

%     figure
%     title('Vm By Direction')
%     hold on
%     plot(outputs.VmByDirection', 'LineWidth', 2)
%     plot(zeros(1, size(params.ExcitationConductances, 2)) + params.SpikeThreshold, '--r', 'LineWidth', 2)
%     xlabel('ms')
%     ylabel('Volts')
%     legend({'0', '45', '90', '135', '180', '225', '270', '315'})
%     
    %% Analyze the results
    
    %Compute the direction selectivity + total tuning curve area +
    %normalized tuning curve area on the basis of peak membrane potential
    peakVmChangeByDirection = max(outputs.VmByDirection - params.RestingVm, [], 2); %Volts
    peakVmChangeByDirection(peakVmChangeByDirection < 0) = 0;
    outputMetricsVm = ConductanceModel_AnalyzeResults(peakVmChangeByDirection', params);
    
    %Compute the direction selectivity + total tuning curve area + normalized tuning curve area
    %on the basis of time spent over threshold (analogous to spikes)
    timeOverThresholdByDirection = sum(outputs.VmByDirection > params.SpikeThreshold, 2)*timestep; %seconds
    %outputMetricsSpikes = ConductanceModel_AnalyzeResults(timeOverThresholdByDirection', params);
    outputMetricsSpikes = ConductanceModel_AnalyzeResults(outputs.SpikesByDirection, params);
    
    DSIs(1, i) = outputMetricsVm.DSI;
    Areas(1, i) = outputMetricsVm.TotalArea;
    NormalizedAreas(1, i) = outputMetricsVm.NormalizedArea;
    VonMisesMeans(1, i) = outputMetricsVm.VonMisesMean;
    VonMisesWidths(1, i) = outputMetricsVm.VonMisesWidth;
    VonMisesErrors(1, i) = outputMetricsVm.VonMisesError;
    PDs(1, i) = outputMetricsVm.PD;
    PDResponses(1, i) = outputMetricsVm.PDResponse;
    NullResponses(1, i) = outputMetricsVm.NullResponse;
    
    DSIs(2, i) = outputMetricsSpikes.DSI;
    Areas(2, i) = outputMetricsSpikes.TotalArea;
    NormalizedAreas(2, i) = outputMetricsSpikes.NormalizedArea;
    VonMisesMeans(2, i) = outputMetricsSpikes.VonMisesMean;
    VonMisesWidths(2, i) = outputMetricsSpikes.VonMisesWidth;
    VonMisesErrors(2, i) = outputMetricsSpikes.VonMisesError;
    PDs(2, i) = outputMetricsSpikes.PD;
    PDResponses(2, i) = outputMetricsSpikes.PDResponse;
    NullResponses(2, i) = outputMetricsSpikes.NullResponse;
end
close(g)

%% Figure out how the above metrics change on the basis of excitation gain
figure
title('DSI')
hold on
xlabel('Excitation Gain')
ylabel('DSI')
plot(ExcitationGainWeights, DSIs(1, :), 'b', 'LineWidth', 2)
semilogx(ExcitationGainWeights, DSIs(2, :), 'r', 'LineWidth', 2)
legend('Vm', 'Spikes')
ylim([0 1])
hold off

figure
title('Area')
hold on
plot(ExcitationGainWeights, Areas(1, :)./max(Areas(1, :)), 'b', 'LineWidth', 2) %normalize to get rid of unit issues and compare the shapes of the curves directly
plot(ExcitationGainWeights, Areas(2, :)./max(Areas(2, :)), 'r', 'LineWidth', 2)
xlabel('Excitation Gain')
ylabel('Area')
legend('Vm', 'Spikes')
hold off

figure
title('Normalized Area')
hold on
plot(ExcitationGainWeights, NormalizedAreas(1, :), 'b', 'LineWidth', 2)
plot(ExcitationGainWeights, NormalizedAreas(2, :), 'r', 'LineWidth', 2)
xlabel('Excitation Gain')
ylabel('Normalized Area')
legend('Vm', 'Spikes')
ylim([0 1])
hold off


% figure
% title('Von Mises Mean')
% hold on
% plot(ExcitationGainWeights, VonMisesMeans(1, :), 'b', 'LineWidth', 2)
% plot(ExcitationGainWeights, VonMisesMeans(2, :), 'r', 'LineWidth', 2)
% xlabel('Excitation Gain')
% ylabel('Von Mises Mu')
% legend('Vm', 'Spikes')
% hold off


% figure
% title('Von Mises Width')
% hold on
% plot(ExcitationGainWeights, VonMisesWidths(1, :), 'b', 'LineWidth', 2)
% plot(ExcitationGainWeights, VonMisesWidths(2, :), 'r', 'LineWidth', 2)
% xlabel('Excitation Gain')
% ylabel('Von Mises Kappa')
% legend('Vm', 'Spikes')
% hold off

% figure
% title('Von Mises Error')
% hold on
% plot(ExcitationGainWeights, VonMisesErrors(1, :)./Areas(1, :), 'b', 'LineWidth', 2) %normalize to response size
% plot(ExcitationGainWeights, VonMisesErrors(2, :)./Areas(2, :), 'r', 'LineWidth', 2) %normalize to response size
% xlabel('Excitation Gain')
% ylabel('Von Mises Error')
% legend('Vm', 'Spikes')
% hold off

figure
title('PDs')
hold on
plot(ExcitationGainWeights, PDs(1, :), 'b', 'LineWidth', 2)
plot(ExcitationGainWeights, PDs(2, :), 'r', 'LineWidth', 2)
xlabel('Excitation Gain')
ylabel('Preferred Direction')
legend('Vm', 'Spikes')
hold off

figure
title('Preferred and Null Responses')
hold on
plot(ExcitationGainWeights, PDResponses(1, :)./max(PDResponses(1, :)), 'b', 'LineWidth', 2) %normalize for units problem. In final figure just use 2 y axes
plot(ExcitationGainWeights, PDResponses(2, :)./max(PDResponses(2, :)), '-r', 'LineWidth', 2)
plot(ExcitationGainWeights, NullResponses(1, :)./max(PDResponses(1, :)), '--b', 'LineWidth', 2)
plot(ExcitationGainWeights, NullResponses(2, :)./max(PDResponses(2, :)), '--r', 'LineWidth', 2)
xlabel('Excitation Gain')
ylabel('Response')
legend('Vm - PD Response', 'Spikes - PD Response', 'Vm - Null Response', 'Spikes - Null Response')
hold off

