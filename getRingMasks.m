function [regMasks,ringrad] = getRingMasks(nreg,dr,ringwidth,celloutline,actcent,imgsize)
% get ring ROIs in the form of image masks
% --------
% input:
% ----------
% nreg = number of ring regions
% dr = radial spacing of regions
% ringwidth = width of rings (set ringwidth=dr to get non-overlapping
% rings)
% celloutline = nx2 coordinates of polygon outlining cell
% imgsize = size of image (rows by columns)
% ------
% output:
% regMasks = list of image masks with individual rings
% -------

cellBW = poly2mask(celloutline(:,1),celloutline(:,2),imgsize(1),imgsize(2));
ringrad = zeros(nreg,1);
for rc = 1:nreg    
    outerrad = ringwidth*3/2 + (rc-1)*dr;
    innerrad = ringwidth/2 + (rc-1)*dr;
    ringrad(rc) = (outerrad + innerrad)/2;
    outerROI = drawcircle(gca,'Center',actcent,'Radius',outerrad,'Visible','off');    
    innerROI = drawcircle(gca,'Center',actcent,'Radius',innerrad,'Visible','off');
    regmask = createMask(outerROI);
    regmaskinner = createMask(innerROI);
    if (rc>1)        
        regmask(logical(regmaskinner)) = 0;
    end
    regMasks(:,:,rc) = regmask.*cellBW;    
end

end