function [newimg,regmask] = removeBrightestRegion(img)
% remove a very bright continuous region in an image
% used to remove overly bright scale bars that make everything else
% difficult to see

% binarize, keeping everything within 99% of the maximal brightness
regmask = imbinarize(img,max(img(:))*0.99);
regmask = bwareafilt(regmask,1); % keep largest connected component
% remove bright region from original image
newimg = img.*(1-regmask);

end