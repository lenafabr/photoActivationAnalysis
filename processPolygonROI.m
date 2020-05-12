function ROIinfo = processPolygonROI(ROI)
% get a structure containing info about an ROI defined by a polygonal
% boundary over an image
% bound is an nx2 set of coordinates

ROIinfo = struct();
ROIinfo.ROI = ROI;
ROIinfo.bound = ROI.Position;     
ROIinfo.mask = createMask(ROI);

% get center of mass and approx radius (treating area as a circle)
polyinfo = polygeom(ROIinfo.bound(:,1),ROIinfo.bound(:,2));
ROIinfo.cent = polyinfo(2:3); % center of mass
% effective radius
ROIinfo.rad = sqrt(polyinfo(1)/pi);

end