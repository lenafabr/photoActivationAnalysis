function [regROIs,cellROI,regdist] = getRandomROIs(img,nreg,regrad,celloutline,actcent,maxdist,mindist)
% generate randomly placed circular ROIs on an image
% the ROIs will be entirely within the polygon defined by celloutline
% returns a list of ROI objects whose position can be manually adjusted
% actcent is the center of the activated region, used to sort the ROIs by
% distance
% maxdist is the maximum allowed distance from actcent
% returns:
% regROIs = adjustable ROI objects for the random regions
% cellROI = adjustable ROI object for the cell
% regdist = distances of each of the regions from the activation center

if (~exist('mindist'))
    mindist = 0;
end

% make mask for the cell ROI
cellBW = poly2mask(celloutline(:,1),celloutline(:,2),size(img,1),size(img,2));
% erode the mask by one region radius
SE = strel('disk',regrad,4)
cellBWer = imerode(cellBW,SE);
%imshowpair(imgs(:,:,end),cellBWer)

% cropped region within which to generate ROIs
%for dc = 1:2
%    regmin(dc) = min(celloutline(:,dc)); regmax(dc) = max(celloutline(:,dc));
%end

% circular region in which to generate ROIs
diffs = celloutline - actcent;
dists = sqrt(sum(diffs.^2,2));
maxR = max(dists);
% maximum allowed distance from the activated region
maxR = min([maxdist,maxR]);
minR = mindist;

% generate random regions within eroded cell mask
% uses rejected sampling: not particularly efficient but good enough
ndone = 0; regcent = [];
while ndone < nreg
    %% generate a bunch of regions
    
    % randomly select center position in circular coords
    newregdist = sqrt(rand(nreg-ndone,1)*(maxR^2 - minR^2) + minR^2);
    newregth = rand(nreg-ndone,1)*2*pi;
    newreg = [cos(newregth),sin(newregth)].*newregdist + actcent;
    

    % get rid of any out of bounds of the image itself
    keepreg = (newreg(:,1)>=1 & newreg(:,1)<=size(img,2) & newreg(:,2)>=1 & newreg(:,2)<=size(img,1));
    newreg = newreg(keepreg,:);
    
    % check which ones fall within the eroded cell mask
    keepreg = cellBWer(sub2ind(size(img),round(newreg(:,2)),round(newreg(:,1))));
    
    % accumulate acceptable regions
    regcent = [regcent; newreg(keepreg,:)];    
    ndone = size(regcent,1);
end

%% sort ROIs by distance from activated region
diffact = regcent-actcent;
regdist = sqrt(sum(diffact.^2,2));
[regdist,sortind] = sort(regdist);

regcent = regcent(sortind,:);

%% display the regions
imshow(img,[])
cellROI = drawpolygon(gca,'Position',celloutline,'Color','w');
hold all
cmap = jet(nreg);
for rc = 1:nreg
    regROIs(rc) = drawcircle('Center',regcent(rc,:),'Radius',regrad,'Label',sprintf('%d',rc),'Color',cmap(rc,:));
end
hold off



end