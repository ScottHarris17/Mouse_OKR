function outputMetrics = ConductanceModel_AnalyzeResults(responseByDirection, params)

orientations = params.orientations; %i.e. stimulus directions

%DSI
[DSI, Vector] = calculateDSI(orientations, responseByDirection);

PD = Vector(1);
%PD/Null Direction Responses
PDResponse = interp1([orientations, 360], [responseByDirection, responseByDirection(1)], PD);
NullDirection = mod(PD + 180, 360);
NullResponse = interp1([orientations, 360], [responseByDirection, responseByDirection(1)], NullDirection);

%Tuning Curve Area
TotalArea = trapz([orientations 360], [responseByDirection, responseByDirection(1)])/360;

%Normalized Tuning Curve Area
NormalizedCurve = responseByDirection./PDResponse;
NormalizedArea = trapz([orientations 360], [NormalizedCurve, NormalizedCurve(1)])/360;

%Von Mises Fit
[fitParams, error] = vonMisesFit(orientations, responseByDirection, deg2rad(90));

%fill the output struct
outputMetrics.DSI = DSI;
outputMetrics.PD = PD;
outputMetrics.PDResponse = PDResponse;
outputMetrics.NullResponse = NullResponse;
outputMetrics.TotalArea = TotalArea;
outputMetrics.NormalizedArea = NormalizedArea;
outputMetrics.VonMisesMean = fitParams(1);
outputMetrics.VonMisesWidth = fitParams(2);
outputMetrics.VonMisesError = error;
end