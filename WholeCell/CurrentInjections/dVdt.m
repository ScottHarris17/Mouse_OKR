ups = grabFromFilter('UP Cells');
downs = grabFromFilter('DOWN Cells');


cutoffFname = 'cutoffValues_lessHyperpolarized.mat'; %choose which cutoff value file to use here

cutoffFilePath = ['C:\Users\mrsco\Box\Completed_Analysis\WholeCell\CurrentInjections\' cutoffFname];
load([cutoffFilePath]);

%choose which protocols to look at
priorityOrder_1 = {'Moving_Bar'}; %first field name shosuld be one of these. Priority is in given order
priorityOrder_2 = {'_SP10'}; %child field name should be one of these. Priority is in given order
mustHave_1 = [];
mustHave_2 = {'CurrentClamp', 'PotassiumSpikes', '_ManyVms'};
cantHave_1 = [];
cantHave_2 = {'Cesium'}; %avoid anything with these tags
rigMandates = []; %no rig mandates

upDvDts_hp = [];
upDvDts_dp = [];
for i = 1:size(ups)
    struct_i = ups{i, 2};
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    
    if pass
        continue
    end
    
    Vms = loc.Analysis_Results.restingMembranePotentialByEpoch;
    cname = struct_i.cellID;
    
    %determine whether cutoff values have been set for this cell yet
    found = 0;
    for j = 1:size(cutoffValues, 1)
        if strcmp(cname, cutoffValues{j, 1})
            found = 1;
            cutoffs = cutoffValues{j, 2};
            break
        end
    end

    hpEpochs = Vms>= cutoffs(1) & Vms <= cutoffs(2); %get hyperpolarazation epochs
    dpEpochs = Vms >= cutoffs(3) & Vms <= cutoffs(4); %get depolarization epochs
    

    hpSpikeTimes = loc.Analysis_Results.allSpikeTimes(hpEpochs);
    dpSpikeTimes = loc.Analysis_Results.allSpikeTimes(dpEpochs);
    
    CE = getClarinetExport(struct_i);
    epochNumsHp = loc.Analysis_Results.EpochNumbers(hpEpochs);
    maxDvDts_hp =[];
    for j = 1:numel(epochNumsHp)
        spikes = hpSpikeTimes{j};
        if numel(spikes) < 1
            continue %skip the epoch if there are no spikes
        end
        epoch_j = CE(epochNumsHp(j)).epoch;
        maxDvDts_hp(end+1) = max(diff(epoch_j))*10000/1000; %volts/s
    end
    upDvDts_hp(end+1) = mean(maxDvDts_hp);
    
    epochNumsDp = loc.Analysis_Results.EpochNumbers(dpEpochs);
    maxDvDts_dp =[];
    for j = 1:numel(epochNumsDp)
        spikes = dpSpikeTimes{j};
        if numel(spikes) < 1
            continue %skip the epoch if there are no spikes
        end
        epoch_j = CE(epochNumsDp(j)).epoch;
        maxDvDts_dp(end+1) = max(diff(epoch_j))*10000/1000; %volts/s
    end
    upDvDts_dp(end+1) = mean(maxDvDts_dp);
end

downDvDts_hp = [];
downDvDts_dp = [];
for i = 1:size(downs)
    struct_i = downs{i, 2};
    [loc, pass] = getLoc(struct_i, priorityOrder_1, mustHave_1, cantHave_1, priorityOrder_2, mustHave_2, cantHave_2, rigMandates);
    
    if pass
        continue
    end
    
    Vms = loc.Analysis_Results.restingMembranePotentialByEpoch;
    cname = struct_i.cellID;
    
    %determine whether cutoff values have been set for this cell yet
    found = 0;
    for j = 1:size(cutoffValues, 1)
        if strcmp(cname, cutoffValues{j, 1})
            found = 1;
            cutoffs = cutoffValues{j, 2};
            break
        end
    end

    hpEpochs = Vms>= cutoffs(1) & Vms <= cutoffs(2); %get hyperpolarazation epochs
    dpEpochs = Vms >= cutoffs(3) & Vms <= cutoffs(4); %get depolarization epochs
    
    hpSpikeTimes = loc.Analysis_Results.allSpikeTimes(hpEpochs);
    dpSpikeTimes = loc.Analysis_Results.allSpikeTimes(dpEpochs);
    
    CE = getClarinetExport(struct_i);
    epochNumsHp = loc.Analysis_Results.EpochNumbers(hpEpochs);
    maxDvDts_hp =[];
    for j = 1:numel(epochNumsHp)
        spikes = hpSpikeTimes{j};
        if numel(spikes) < 1
            continue %skip the epoch if there are no spikes
        end
        epoch_j = CE(epochNumsHp(j)).epoch;
        maxDvDts_hp(end+1) = max(diff(epoch_j))*10000/1000; %volts/s
    end
    downDvDts_hp(end+1) = mean(maxDvDts_hp);
    
    epochNumsDp = loc.Analysis_Results.EpochNumbers(dpEpochs);
    maxDvDts_dp =[];
    for j = 1:numel(epochNumsDp)
        spikes = dpSpikeTimes{j};
        if numel(spikes) < 1
            continue %skip the epoch if there are no spikes
        end
        epoch_j = CE(epochNumsDp(j)).epoch;
        maxDvDts_dp(end+1) = max(diff(epoch_j))*10000/1000; %volts/s
    end
    downDvDts_dp(end+1) = mean(maxDvDts_dp);
end

figure
title('Max dV/dt - depol. vs hyperpol.')
hold on
scatter(upDvDts_dp, upDvDts_hp, 'k', 'filled', 'MarkerFaceAlpha', 0.6)
scatter(downDvDts_dp, downDvDts_hp, 'm', 'filled', 'MarkerFaceAlpha', 0.6)
plot([0 200], [0 200], '--r')
xlabel('Depolarized')
ylabel('Hyperpolarized')
p = signrank([upDvDts_dp, downDvDts_dp], [upDvDts_hp, downDvDts_hp])
signrank(upDvDts_dp, upDvDts_hp)
signrank(downDvDts_dp, downDvDts_hp)
