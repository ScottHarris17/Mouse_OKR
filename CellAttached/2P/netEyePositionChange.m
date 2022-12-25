function netChange = netEyePositionChange(upPSTHs, downPSTHs)
%computes net change in eye position given oDSGC PSTHs to oscillating
%grating
if size(upPSTHs, 1)>1
    medUp = mean(upPSTHs);
else
    medUp = upPSTHs;
end

if size(downPSTHs, 1)>1
    medDown = mean(downPSTHs);
else
    medDown = downPSTHs;
end

velocity = medDown-medUp;
position = cumtrapz(velocity);
netChange = position(end) - position(1);
end
