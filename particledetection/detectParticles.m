function [goodpts,imgclean, pts] = detectParticles(imgorig,w,varargin)
% -------------
%% detect bright particle locations in an image
% w is expected particle radius
% img should be loaded in with imread
% also returns cleaned up image
% maxecc: set a maximum eccentricity cutoff
% intrat: minimum intensity per pixel ratio
% massmin and massmax are minimal and maximal feature masses 
% Imax: set a minimal peak brightness cutoff in the cleaned up image (top
% percent of cleaned image)
% given as top percentile to keep
% dodisplay: display image
% optional: if fluorescent is set on then:
    % a) image is not flipped
    % b) bandpass filter rather than gaussian filter is used
% outputs all retained points with the following info (f=goodpts)
% 		f(:,1):	the x centroid positions, in pixels.
% 		f(:,2): the y centroid positions, in pixels. 
% 		f(:,3): integrated brightness of the features. ("mass")
% 		f(:,4): the square of the radius of gyration of the features.
% 		    (second moment of the "mass" distribution, where mass=intensity)
% 		f(:,5): eccentricity, which should be zero for circularly symmetric features and 
%                   order one for very elongated images

Imax = 0.3;
massmin = 0;
massmax = inf;
maxecc = 0.5;
intrat = 0;
dodisplay = 0;
flipimage = 0;
filtertype = 'none';

lambda = 1; 

if (~isempty(varargin))
    for vc = 1:2:length(varargin)
        switch varargin{vc}
            case('Imax')
                Imax = varargin{vc+1};
            case('massmin')
                massmin = varargin{vc+1};
            case('massmax')
                massmax = varargin{vc+1};
            case('maxecc')
                maxecc = varargin{vc+1};
            case('intrat')
                intrat = varargin{vc+1};
            case('dodisplay')
                dodisplay = varargin{vc+1};
            case('flipimage')
                flipimage = varargin{vc+1};
            case('filtertype')
                filtertype = varargin{vc+1};
        end
    end
end
    
% clean up image with some band pass or low pass filter
if (flipimage)    
    img = max(max(imgorig))-imgorig;
else
    img = imgorig;
end

% clean up image with some band pass or low pass filter
if (strcmp(filtertype,'gauss'))
    % clean up with gaussian filter
    G = fspecial('gaussian',[w w],w/2);
    imgclean = imfilter(img,G,'same');    
elseif (strcmp(filtertype,'bpass'))
    % clean up with bandpass filter   
    imgclean = bpass(img,lambda,w);    
elseif (strcmp(filtertype,'none'))
    imgclean = imgorig;
else
    error('Unknown filtertype.', filtertype)        
end

% find features starting with clean image
pts = feature2D(imgclean,lambda,w,massmin,Imax,2,'none');

if (dodisplay==1)
    imshow(imgorig,[]);
elseif (dodisplay==2)
    imshow(img,[]);
elseif (dodisplay==3)
    imshow(imgclean,[]);
end

if (isempty(pts))
    goodpts = [];
    return
end

goodind = find(pts(:,5)<maxecc & pts(:,3)./(pi*pts(:,4))>intrat & pts(:,3)<massmax);
goodpts = pts(goodind,:);



if (dodisplay>0)
    hold all
    plot(goodpts(:,1),goodpts(:,2),'b.','MarkerSize',5)
    hold off
end


