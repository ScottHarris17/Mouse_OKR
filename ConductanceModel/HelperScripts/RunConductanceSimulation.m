function outputs = RunConductanceSimulation(params, timestep)
%This function runs the simulation using the differential equation for
%membrane potential and integrating with Euler's method

%differential equation for dV/dt
dvdt = @(Vm, Gex, Ex, Gin, Ein, Grest, Erest, C)...
     (Gex*(Ex-Vm) + Gin*(Ein-Vm) + Grest*(Erest-Vm))/C;

% %example of how to test the effect of a current injection on the model cell
% exStep = zeros(1, 5500);
% exStep(1000:3000) = -10;
% Gex = 1e-12*exStep./(params.ExcitationHold - params.ExcitationReversal);

%Grab the Parameters that are needed from the input structure
orientations = params.orientations; %i.e. stimulus directions
Ex = params.ExcitationReversal;
Ein = params.InhibitionReversal;
Grest = 1/params.InputResistance;
Erest = params.RestingVm;
C = params.MembraneCapacitance;

%initialize
VmByDirection = zeros(numel(orientations), size(params.ExcitationConductances, 2));
VmPlusSpikes = zeros(numel(orientations), size(params.ExcitationConductances, 2));
spikesByDirection = zeros(1,numel(orientations));

for i = 1:numel(orientations)
    Vm = params.RestingVm;
    Gex = params.ExcitationConductances(i, :)*params.ExcitationGain(i) * params.ExcitationGainWeight;
    Gin = params.InhibitionConductances(i, :)*params.InhibitionGain(i) * params.InhibitionGainWeight;
    
    %use Euler's method to integrate the differential equation at each
    %subsequent time step
    for j = 1:size(VmByDirection, 2)
        VmByDirection(i, j) = Vm;
        Gex_now = Gex(j);
        Gin_now = Gin(j);
        deltaV = dvdt(Vm, Gex_now, Ex, Gin_now, Ein, Grest, Erest, C)*timestep;
        Vm = Vm + deltaV;
    end
    
    %repeat the simulation,  but now calculate spikes
    j = 1;
    while j <= size(VmPlusSpikes, 2)
        if Vm < params.SpikeThreshold % no spike
            VmPlusSpikes(i, j) = Vm;
            Gex_now = Gex(j);
            Gin_now = Gin(j);
            deltaV = dvdt(Vm, Gex_now, Ex, Gin_now, Ein, Grest, Erest, C)*timestep;
            Vm = Vm + deltaV;
            j = j + 1;
        else
            %if there's a spike note it, skip 3 ms and return to resting Vm
            spikesByDirection(i) = spikesByDirection(i) + 1;
            j = j+3;
            if j > size(VmPlusSpikes, 2)
                break
            end
            Vm = params.RestingVm; %reset membrane potential after each spike
        end       
    end
end

outputs.VmByDirection = VmByDirection;
outputs.VmPlusSpikes = VmPlusSpikes;
outputs.SpikesByDirection = spikesByDirection;
end