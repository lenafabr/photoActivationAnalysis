%% This is an example script for loading data from a single movie for 
% a particular cell.
% The results are saved to a cell object for future use
% saves one image frame, ROI definitions, and time traces only
% so should not generate very large data files

% To process many cells at once more efficiently, use
% exampleManyCellProcessing.m


%% establish input/output file names
cellname = 'Cecile_WT_FRAP001'; % short name associated with his movie or cell

% directory where raw images are stored
dirname = '/data/proj/ERtransport/ER_paGFP_Cecile/ER_paGFP_Cecile/COS7-WT-paGFP-KDEL-ERmCherry-040822/FRAP 001/';

% where to save the final .mat file with this cell object
savefile = sprintf('../celldata/example_%s.mat',cellname);

% get file names for images to load in:
% 1) photoactivated channel, prior to photoactivation
prefile = 'FRAP 001_FRAP Pre Series01_ch02.tif';
% 2) photoactivated signal channel, during photoactivation
PAfile = 'FRAP 001_FRAP Bleach Series03_ch02.tif';
% 3) photoactivated region only, during photoactivation
PAregfile = 'FRAP 001_FRAP Bleach Series03_ch00.tif';
% 4) general ER luminal marker, during photoactivation
ERfile = 'FRAP 001_FRAP Bleach Series03_ch03.tif';
% 5) metadatafile exported from ImageJ (used to get frame rate)
% (Note: Lena exports this from ImageJ when reading in .lif files
% it should be possible to parse the frame rate out of the .xml metadata 
% files as well);
% currently ignore
metafile='';

%% create a cell object and load data into it
% Note: this does NOT load the actual raw image files (except for one frame
% saved as CL.ERimg saved for quick visualization
CL = CellObjPA(cellname);
if (isempty(metafile))
    CL.loadCellData(dirname,PAfile,prefile,PAregfile,ERfile);
else
    CL.loadCellData(dirname,PAfile,prefile,PAregfile,ERfile,metafile);
end

% FOR CECILE'S DATA only: 
% blank out the scale bar as it makes everything else difficult to see
[CL.ERimg,sclbarmask] = removeBrightestRegion(CL.ERimg);
imshow(CL.ERimg,[])
title(sprintf('%s',CL.Name),'Interpreter','none')


%% Conversion factors to real length and time units
% Ought to get this from metadata, but for now 
% SET THIS MANUALLY:
sclbarlen = 10; % what is the scale bar length in um? 
CL.dt = 1; % seconds per frame

% get the px per um from the scalebar
barstats = regionprops(sclbarmask);
% WARNING: this might be off by 1 pixel depending on how the bar was drawn
CL.pxperum = barstats.BoundingBox(3)/sclbarlen;

%% manually identify rectangular activation region
CL.getActROI(struct('dodisplay',1,'manualrect',1));
    
%% manually identify cell and nucleus outlines, based on the general ER region

% for comparison also show photoactivated signal at end of photoactivation period
img = imread([CL.DirName CL.PAfile], CL.NFrame-CL.startPA-1);
img = imadjust(img,[0,0.3],[0,1],0.5);
figure(2)
imshow(img)
title(sprintf('%s: photoactivated signal',CL.Name),'Interpreter','none')
%
figure(1)
CL.setCellNucROI(struct('showPA',0));

% outer boundary of region to be analyzed set to the extracted cell
% boundary
CL.cellROI = CL.fullcellROI;

%% optional: show and adjust cell and nucleus segmentation as needed
CL.showCellROI(true);

%% save object to file
save(savefile,'CL')

%% establish wedge ROIs for tracking signal at different distances from photoactivation center
% this step can be a somewhat slow

dR = 1; % radius separations (in um)
Rwidth = 2; % thickness of rings (in um);
maxR = 15; % maximal outer radius (in um)    

optionsROIs = struct();
optionsROIs.dodisplay = 2; % draw all regions
optionsROIs.minR = 2; % minimal inner radius (in um)  
optionsROIs.arclen = 2; % arc length in um for each wedge
optionsROIs.arcshift = 1; % arc length shift in um between neighbor wedges
% minimum length (um) of ray to cell boundary in order to keep wedge;
% this parameter selects wedges in a bulk region of the cell
% not squeezed between photoactivation center and nearby cell boundary
% not on the other side of the nucleus from the PA center
optionsROIs.rayR = 5; 
% how far (in um) to erode the cell bulk segmented region (to avoid regions
% right next to cell boundary)
optionsROIs.erodemask = 1;
% keep only wedges with at least this fraction of their area included in
% the segmented cell region
optionsROIs.minareafrac = 0;
% keep only wedges with this minimal area (in um^2) included in the
% segmented cell region
optionsROIs.minarea = 2;
% also get the whole-ring ROIs as well as the wedges
optionsROIs.getRingROIs = true;

% Rvals are the radii at the center of each ring band
[Rvals,whichrad] = CL.getWedgeROIs(maxR,dR,Rwidth,optionsROIs);


%% get time-traces over each region
% get non-photoactivated luminal marker traces?
getnonPATrace = false;
% imgs saves all raw images
% regionTraces saves traces of fluorescence intensity per area for each
% wedge or ring ROI
[regionTraces,imgs] = CL.getROItraces(getnonPATrace);  
 

%% Plot example ROIs and traces to make sure everything makes sense
% This is just for visualization
% angle (in degrees, btwn -180 and 180) around the circle, indicating roughly
% which wedge to show for each ring
% 0 degrees points east. Angle increases clockwise.
pickangle = -20; 

CL.plotExampleWedgeROIs(pickangle)

%% Fit signal traces to a smoothed double-exponential function
optionsFit = struct();
% use a double exponential fit, with max value set by the max of the
% photoactivated region.
optionsFit.fittype = '2expfixlim';
% get the time to half-max
optionsFit.cutfrac = 0.5;
% avoid very noisy traces:
% 1) throw out all ROIs where the signal in the last few frames 
% is below this many standard deviations of the pre-activation signal
optionsFit.minsignalchange = 1.96;
% 2) throw out all ROIs where the signal in the last few frames
% has a standard deviation over mean above this value
optionsFit.maxendfrano = 1;

halftimes = CL.getHalfTimes(optionsFit);

%% get median signal for each distance and plot to see scaling
whichrad = [CL.ROIs.whichrad];
iswedge = cellfun(@(x) contains(x,'wedge'),{CL.ROIs.type});
medhalftime = zeros(1,length(Rvals));
for rc = 1:length(Rvals)
    wedgeind = find(iswedge & whichrad==rc);
    
    htimes = [CL.ROIs(wedgeind).halftime];
    
    medhalftime(rc) = nanmedian(htimes);
end

figure
loglog(Rvals,medhalftime,'.-')
hold all
loglog(Rvals,Rvals.^2,'--')
loglog(Rvals,10*Rvals,'--')
legend('results','R^2','10R')
xlabel('distance (um)')
ylabel('signal rise time (sec)')
set(gca,'FontSize',14)
title(CL.Name,'Interpreter','none')

%% save object to file
save(savefile,'CL','Rvals','optionsROIs','optionsFit')
