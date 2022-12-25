function [centers, returnImage] = findCenters_MosaicAnalysis(img, threshold, localRadius, blankArea, blurSize, minimumArea, pixelsPerMM)

%% Convert to pixels
localRadius = round((localRadius/1000)*pixelsPerMM);
minimumArea = round((minimumArea/1000)*pixelsPerMM);
blankArea = round((blankArea/1000)*pixelsPerMM);

close all
bw = rgb2gray(img);

%compute local averages
SE=ones(localRadius,localRadius);SE=SE/(sum(SE(:))-blankArea);
SE(round(size(SE,1)/2-sqrt(blankArea)/2):round(size(SE,1)/2+sqrt(blankArea)/2),...
   round(size(SE,2)/2-sqrt(blankArea)/2):round(size(SE,2)/2+sqrt(blankArea)/2))...
   = 0; %blank out the center of the kernel
localAverages = uint8(conv2(bw,SE,'same')); %convolve with a kernel to get local average

%blur the image to get rid of granularity
blur = imgaussfilt(bw, blurSize);

%threshold based on standard deviation (computed once for the image, not on
%a local basis)
er = std(double(blur(:)));
t = blur>(localAverages+threshold*er);

figure(1);
imshow(t)

%user can remove cells (or parts of cells) if needed
while 1
    resp = questdlg('Edit?');
    if ~strcmp(resp, 'Yes')
        break
    end
    
    rect = drawrectangle;
    p = round(rect.Position);
    
    t(p(2):p(2)+p(4), p(1):p(1)+p(3)) = 0;
    figure(1);
    imshow(t)
end

boundaries = bwboundaries(t); %group contiguous areas of white
centers = []; %write down the center coordinates
for i = 1:numel(boundaries)
    b_i = boundaries{i};
    
    xCords = b_i(:, 1);
    yCords = b_i(:, 2);
    area = polyarea(xCords, yCords);
    if area < minimumArea %throw out object if too small
        continue
    end
    centers = [centers;mean(b_i)];
    
    figure(1)
    hold on
    plot(b_i(:, 2), b_i(:, 1));
end
hold off

figure(2)
imshow(img)
hold on
scatter(centers(:, 2), centers(:, 1), 'ow')
hold off
%user can add cells if needed
while 1
    resp = questdlg('Edit Centers?', 'Editor', 'Add', 'Remove', 'Done', 'Add');
    if strcmp(resp, 'Done')
        break
    end
    
    if strcmp(resp, 'Add')
        newPoint = drawcrosshair();
        p = round(newPoint.Position);

        centers = [centers; p(2), p(1)];
    
    elseif strcmp(resp, 'Remove')
        rect = drawrectangle();
        p = round(rect.Position);
        removing = [];
        for i = 1:size(centers, 1)
            if centers(i, 1) > p(2) && centers(i, 1) < p(2) + p(4) &&...
                    centers(i, 2) > p(1) && centers(i, 2) < p(1) + p(3)
                removing = [removing, i];
            end
        end
        centers(removing, :) = [];
    end
    returnImage = figure(2);
    imshow(img)
    hold on
    scatter(centers(:, 2), centers(:, 1), 'ow')
    hold off
end
returnImage = figure(2);
imshow(img)
hold on
scatter(centers(:, 2), centers(:, 1), 'ow')
hold off     
end