function [indices, saccadeVelocities] = findSaccades_old(trace, velocityThreshold, FrameRate)
%Calculates start/stop indices and velocity of each saccade
%inputs:
%     - trace: vector of eye positions in degrees
%     - velocityThreshold: minimum deg/sec eye movements to classify as saccades
%     - FrameRate: frame rate of the video recording
% outputs:
%     - indices: nx2 matrix where n is the number of saccades detected.
%             The first column is the start index of each saccade and the second
%             column is the end index
%     - saccadeVelocities: velocity (direction and speed) of the first frame of each saccade     

minInterval = 1200; %ms - minimum time between two different saccades
minNumberOfFrames = ceil((minInterval/1000)/(1/FrameRate));

eye_velocity = diff(trace);

indices_prelim = find(abs(eye_velocity) > velocityThreshold/FrameRate);

%remove book ended indices
indices_prelim = indices_prelim(indices_prelim > 10 & indices_prelim < (numel(trace) - 10));

%Check edge case where there are no saccades
if numel(indices_prelim) == 0
    indices = [];
    saccadeVelocities = [];
    return
end

%figure out start times of each saccade
startIndices = [indices_prelim(1)];
saccadeVelocities = [eye_velocity(indices_prelim(1))];
for i = 1:numel(indices_prelim)-1
    saccadeInterval = indices_prelim(i+1) - indices_prelim(i);
    if saccadeInterval < minNumberOfFrames
        continue
    else
        startIndices(end +1, 1) = indices_prelim(i+1);
        saccadeVelocities(end + 1, 1) = eye_velocity(indices_prelim(i + 1));
    end
end

%figure out end times of each saccade
endIndices = [];
%add the end frame of every saccade
for i = 2:numel(indices_prelim)
    saccadeInterval = indices_prelim(i) - indices_prelim(i-1);
    if saccadeInterval < minNumberOfFrames
        continue
    else
        endIndices(end + 1, 1) = indices_prelim(i-1)+1; %here you're marking the fram after every saccade, so add 1
    end
end
endIndices(end + 1, 1) = indices_prelim(end) + 1;

indices = [startIndices, endIndices];
end
