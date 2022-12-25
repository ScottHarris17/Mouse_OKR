function output = VmCurrentInjection_Helper(traces, directions, varargin)
%called by VmCurrentInjection.m.
% Inputs:
%   - traces = nxm matrix. n = observation number, m = trace number
%   - directions = stimulus direction for each trace
% Outputs:
%   - structure with the DSI, Area, and normalized area of the tuning curve
%   generated from the traces and their corresponding stimulus directions.
%   Contains metrics calculated using area under the Vm trace (AUCs) and
%   the peak Vm response (MaxResponse)


%optional input parameters:
%     - Offset = a single number to add (or subtract) from each max Vm value (e.g. for ciruclarization)
%     - threshold = a single number below which max Vm values are thresholded back to 0

defaultOffset = 0;
defaultThreshold = 0;

p = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x);
addParameter(p,'offset',defaultOffset, validScalar);
addParameter(p, 'threshold', defaultThreshold, validScalar);
parse(p, varargin{:});
offset = p.Results.offset;
threshold = p.Results.threshold;

   
allDirections = [0:45:315];
restingVmByEpoch = zeros(1, numel(allDirections));

allTracesByDirection = cell(1, numel(allDirections));
for i = 1:numel(directions)
    trace_i = traces(:, i)';
    restingVm = mean(trace_i(1:1000));
    adjustedTrace_i = trace_i - restingVm;
    direction_i = directions(i);
    thisMat = allTracesByDirection{allDirections == direction_i};
    newMat = [thisMat; adjustedTrace_i];
    allTracesByDirection{allDirections == direction_i} = newMat;
    restingVmByEpoch(i) = restingVm;
end
meanByDirection = zeros(numel(allDirections), size(traces, 1));
SEMByDirection = zeros(numel(allDirections), size(traces, 1));
for i = 1:numel(allDirections)
    traces_i = allTracesByDirection{i};
    meanByDirection(i, :) = mean(traces_i);
    SEMByDirection(i, :) = std(traces_i)/sqrt(size(traces_i, 1));
end

meanByDirection(isnan(meanByDirection)) = 0; %sometimes there will be trailing NaNs due to interpolation, change to 0's
meanByDirectionSmoothed = movmean(meanByDirection, 10, 2);
AUCs = trapz(meanByDirectionSmoothed, 2)';
MaxResponses = max(meanByDirectionSmoothed, [], 2)';


%% Analysis for both AUCs and max responses:
[DSI_AUCs, Vector_AUCs] = calculateDSI(allDirections, AUCs);

%PD/Null Direction Responses
PD_AUCs = mod(Vector_AUCs(1), 360);
PDResponse_AUCs = interp1([allDirections, 360], [AUCs, AUCs(1)], PD_AUCs);
NullDirection_AUCs = mod(PD_AUCs + 180, 360);
NullResponse_AUCs = interp1([allDirections, 360], [AUCs, AUCs(1)], NullDirection_AUCs);

%Tuning Curve Area
TotalArea_AUCs = trapz([allDirections 360], [AUCs, AUCs(1)])/360;

%Normalized Tuning Curve Area
NormalizedCurve_AUCs = AUCs./PDResponse_AUCs;
NormalizedArea_AUCs = trapz([allDirections 360], [NormalizedCurve_AUCs, NormalizedCurve_AUCs(1)])/360;

output.DSI_AUCs = DSI_AUCs;
output.PD_AUCs = PD_AUCs;
output.PDResponse_AUCs = PDResponse_AUCs;
output.NullResponse_AUCs = NullResponse_AUCs;
output.TotalArea_AUCs = TotalArea_AUCs;
output.NormalizedArea_AUCs = NormalizedArea_AUCs;


%repeat for max responses
MaxResponses = MaxResponses + offset; %offset is usually 0 but see varargin above
MaxResponses(MaxResponses < threshold) = 0; %threshold to 0


[DSI_MaxResponses, Vector_MaxResponses] = calculateDSI(allDirections, MaxResponses);

PD_MaxResponses = mod(Vector_MaxResponses(1), 360);
PDResponse_MaxResponses = interp1([allDirections, 360], [MaxResponses, MaxResponses(1)], PD_MaxResponses);
NullDirection_MaxResponses = mod(PD_MaxResponses + 180, 360);
NullResponse_MaxResponses = interp1([allDirections, 360], [MaxResponses, MaxResponses(1)], NullDirection_MaxResponses);


%Tuning Curve Area
TotalArea_MaxResponses = trapz([allDirections 360], [MaxResponses, MaxResponses(1)])/360;

%Normalized Tuning Curve Area
NormalizedCurve_MaxResponses = MaxResponses./PDResponse_MaxResponses;
NormalizedArea_MaxResponses = trapz([allDirections 360], [NormalizedCurve_MaxResponses, NormalizedCurve_MaxResponses(1)])/360;

if sum(MaxResponses) == 0 %in the case that you get no responses, perfectly DS (usually only caused by thresholding)
    DSI_MaxResponses = 1;
    NormalizedArea_MaxResponses = 0;
    TotalArea_MaxResponses = 0;
end

output.MeanRestingMembranePotential = mean(restingVm);
output.DSI_MaxResponses = DSI_MaxResponses;
output.PD_MaxResponses = PD_MaxResponses;
output.PDResponse_MaxResponses = PDResponse_MaxResponses;
output.NullResponse_MaxResponses = NullResponse_MaxResponses;
output.TotalArea_MaxResponses = TotalArea_MaxResponses;
output.NormalizedArea_MaxResponses = NormalizedArea_MaxResponses;
output.MaxResponsesByDirection = MaxResponses;
output.MaxResponseOffset = offset;
output.Threshold = threshold;
output.offset = offset;