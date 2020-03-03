function regionTraces = getRegionTraces(imgs,regROIs,regMasks)
% get brightness over time for each ROI region
% the region ROIs must exist (do not delete the figure!)
% -----
% Input
% -------
% imgs = (height x width x frames) image series
% regROIs = array of ROI objects
% optionally: supply regMasks (masks for each region) and can leave regROIs
% as undefined
% -------
% Output
% ------
% returns a matrix of brightness per unit area in each region at each time
% point
% regionTraces(rc,fc) is for region rc, time frame fc

if (exist('regMasks','var'))
    % use predefined region masks
    allmasks = regMasks;
else
% create a mask for each ROI
for rc = 1:length(regROIs)
    circmask = createMask(regROIs(rc));    
    allmasks(:,:,rc) = circmask;
end
end

% area in each mask (in px^2)
areas = squeeze(sum(sum(allmasks,1),2));

% reshape each image into a single vector
imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
maskshape = reshape(allmasks,size(imgs,1)*size(imgs,2),size(allmasks,3));

% get matrix of total brightness for each region, for each frame
% brightness is scaled by the region area
%regrads = [regROIs.Radius]';
regionTraces = (maskshape'*imgshape)./areas;

end