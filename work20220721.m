% general directory name
upperdirname = '/media/ekoslover/PortableDrive2/pc1/.snapshots/snapshot.0/data/proj/ERtransport/calcium_release/';
gendirname = [upperdirname 'COS7-WT-GCaMP3-cIP3-Flash-210422/'];
%gendirname = [upperdirname 'COS7-ATLKO-gcamp3-cIP3-short-pulse-280422'];

% cell name
cellname = '210422-004-puff';

% specific directory for this cell
celldir= sprintf('%s/',cellname);

% filenames for pre, bleach, post-bleach

% before uncaging
prefile = '';
    
% post-uncaging
fglob = sprintf('%s*Pb1*_ch00.tif',cellname);
tmp = dir([gendirname celldir fglob]);
PAfile = tmp.name;

%% create a cell object and load data    
CL = CellObjPA(cellname);
    
CL.loadCellData([gendirname celldir],PAfile,prefile,PAfile,PAfile);

CL.ERimg = double(CL.ERimg)/max(double(CL.ERimg(:)));

%% load in all images
CL.imgs = loadImages(CL.DirName,CL.PAprefile,CL.PAfile);

%% look at image to make sure it looks ok
imshow(CL.imgs(:,:,1),[0,0.3])


%%
addpath('./roisignalGUI/')

%% test gui
app = roiSignalProcess('CL',CL)



%% --------------- play with spatial peak finding ------------
allmasks = compmask;
imgs = CL.imgs;
% reshape each image into a single vector
imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
maskshape = reshape(allmasks,size(imgs,1)*size(imgs,2),size(allmasks,3));

%% get matrix of total brightness for each region, for each frame
% brightness is scaled by the region area
%regrads = [regROIs.Radius]';
regionTrace = (maskshape'*imgshape);
tvals = 1:length(regionTrace);
%%
plot(tvals,regionTrace)

%% 
newROIs = [ROIs(1)];

newROIs(1).mask = compmask;
newROIs(1).avgsignal = regionTrace;

options = struct();
options.detrendorder = str2num(app.DetrendorderEditField.Value);
options.order = str2num(app.OrderEditField.Value);
options.framelen = str2num(app.WindowEditField.Value);
options.cutscalepeak = str2num(app.PeakthresholdEditField.Value);
options.cutscalecurvature = str2num(app.CurvaturethresholdEditField.Value);


[keepreg, newROIs] = filterForPeaks(newROIs,options);

%% plot signal in entire region
plot(tvals,newROIs.avgsignal)
hold all
puffind = newROIs.puffind;
plot(tvals(puffind),newROIs.avgsignal(puffind),'o')
hold off

%% plot image at peak time
figure
imshowpair(CL.imgs(:,:,puffind(1)),compmask)

%% look for spatial peak (treat as particle finding problem)
img = CL.imgs(:,:,puffind(1)+2);
peaksize = 5; % size of peak to be used in bandpass filter
% dilate mask
%se = offsetstrel('ball',15,15);
se = strel('disk',peaksize*2)
dilmask = imdilate(compmask,se);
imshowpair(img,dilmask)
imgmask = dilmask.*img;
imgmask0 = compmask.*img;
imshow(imgmask)

%% particle finding options
% default options
opt = struct();
% expand segmentation contour by a given distance
% and search only particles within the expanded contour
% if empty then search for all particles within preset cell rectangle
opt.expandcont = [];
% amount of displaying (0 is no figures drawn, >1 does
% intermediate steps)
opt.dodisplay = 2;

% particle finding parameters
% upper limit width (diameter) of features, used in bandpass filter
opt.w = peaksize;%ceil(0.6*pxperum);
% intensity ratio cutoff
% only keep features with a brightness per unit area that's at least this
% fraction of the maximum (with some top bright outliers cut out to
% calculate maximum)
opt.intratscl = 0;
% percentile cutoff for overall feature mass
opt.ppmass = 95;

% maximal eccentricity
opt.maxecc = inf;

% invert image to find dark particles
opt.invert = 0;
% limits for contrast stretching
% percentiles to min out and max out
opt.stretchlim=[0,1];

% brightness scaling for displaying image
opt.showscale = [];

% type of average to use in brightness calculations
opt.bravgtype = 'median';

% crop to integer pixel values
opt.intpixcrop = 1;

partopt = opt;

% adjust image contrast
imgA = imadjust(imgmask,stretchlim(imgmask,opt.stretchlim));

% bandpass filter
imgclean = bpass(imgA,1,opt.w);
maxclean = max(imgclean(:));

% parameters for feature finding
Imax= 0.0025;
massmin = prctile(imgclean(:),opt.ppmass)*(2*opt.w)^2;
massmax = inf;
% allow for a bright outlier
pp = 1-(2*opt.w+1)^2/prod(size(imgclean));
maxclean = prctile(imgclean(:),pp*100);
intrat = opt.intratscl*maxclean;

[goodpts,tmp,pts] = detectParticles(imgclean,opt.w,'Imax',Imax,'massmin',massmin,'massmax',massmax,...
    'maxecc',opt.maxecc,'intrat',intrat);

% keep only those points within the original mask
maskvals = interp2(compmask,goodpts(:,1),goodpts(:,2));
goodpts = goodpts(abs(maskvals-1)<eps,:);

%% plot found peaks
%imshowpair(imgclean,compmask)
imshowpair(img,compmask)
%imshow(img,[])
hold all
plot(goodpts(:,1),goodpts(:,2),'m.','MarkerSize',10)
%[~,topind] = max(goodpts(:,3)) % find biggest peak
goodpts = sortrows(goodpts,3,'descend');
topind = 1
plot(goodpts(topind,1),goodpts(topind,2),'go','MarkerSize',10)
for pc = 1:size(goodpts,1)
    text(goodpts(pc,1),goodpts(pc,2),sprintf('%d',pc),'Color','c')
end
hold off

% -----------------------------------
%% time smoothing of images
% --------------------------------------

% make a gaussian window filter
w = ceil(newROIs(1).puffwidths(1)/2);% half-width window
gw = gausswin(w*2+1,w);
% get weighted time average of ROIs;
inds = puffind(1)-w:puffind(1)+w;
permimg= permute(CL.imgs(:,:,inds),[3,1,2]);
imgtimefilt = squeeze(pagemtimes(gw',permimg));




%% smooth image with averaging filter
Kaverage = filter2(fspecial('average',10),imgtimefilt)/255;
%imshow(Kaverage,[])
imgsmoothed = Kaverage.*dilmask;
imshow(imgsmoothed,[])
%imgsmoothed = bpass(imgmask,1,20);
%imgsmoothed = imgsmoothed.*compmask;

[goodpts,tmp,ptcleans] = detectParticles(imgsmoothed,10,'Imax',Imax,'massmin',massmin,'massmax',massmax,...
    'maxecc',opt.maxecc,'intrat',intrat);
% keep only those points within the original mask
maskvals = interp2(compmask,goodpts(:,1),goodpts(:,2));
goodpts = goodpts(abs(maskvals-1)<eps,:);
goodpts = sortrows(goodpts,3,'descend');

hold all
plot(goodpts(:,1),goodpts(:,2),'m.','MarkerSize',10)
hold off

[~,maxind] = max(imgsmoothed(:));
[y,x] = ind2sub(size(imgsmoothed),maxind);
hold all
plot(x,y,'co')
hold off

goodpts = sortrows(goodpts,3,'descend');
peakcent = goodpts(1,1:2)
%% show final peak location
imshowpair(imgtimefilt,compmask)
%imshow(imgtimefilt,[])
hold all
plot(peakcent(1),peakcent(2),'c.','MarkerSize',10)
hold off

%% make an ROI around the peak
roiwidth = 20;
rect = [peakcent-[roiwidth/2,roiwidth/2] roiwidth roiwidth];
imshow(imgtimefilt,[])
hold all
plot(peakcent(1),peakcent(2),'m.','MarkerSize',10)
peakroi = drawrectangle('Position',rect,'Color','g')
peakroimask = createMask(peakroi);
hold off

%% plot signal for ROI
allmasks = peakroimask;
% reshape each image into a single vector
imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
maskshape = reshape(allmasks,size(imgs,1)*size(imgs,2),size(allmasks,3));

regionTracePeak = (maskshape'*imgshape);
tvals = 1:length(regionTracePeak);

plot(tvals,regionTracePeak)
hold all
plot(tvals(puffind(1)),regionTracePeak(puffind(1)),'ro')
hold off
%% find center of mass for cleaned masked image (not useful)
props = regionprops(compmask, imgsmoothed, 'Centroid', 'WeightedCentroid');
com = props.Centroid;
hold all
plot(com(1),com(2),'y*','LineWidth',2)
