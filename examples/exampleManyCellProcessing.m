% -------
%% Example file for processing data from many experimental images
% To look at or adjust a single cell at a time, look at exampleCellProcessing.m
% --------

% Which file numbers do you want to process?
% by default, for number 1, this will look for files that have FRAP001 in
% the name, etc. But this can be adjusted
processnums = [1,2];

%% establish input/output file names
% pattern for name associated with each cell object
cellnamepattern = 'COS7_WT_FRAP%03d'; 

% directory where raw images are stored
dirname = '/data/proj/ERtransport/210423_COS7_ER-PAGFP_RTN4KO/WT/';

% where to save .mat file with all the cell objects together
savefile = '../celldata/example_WT.mat';

% first part of file name, that all input iamge files have in common
genfilename = '210423_COS7_ER_PAGFP_ER_mCherry_WT_';

% get file names for images to load in:
% these are currently set up as file name globs on a linux system
% might need to adjust for other OS;
% 1) photoactivated channel, prior to photoactivation
prefileglob = 'FRAP_%03d_*_Pre_*C2.tif';
% 2) photoactivated channel, during photoactivation
PAfileglob = 'FRAP_%03d_*_Bleach_*C2.tif';
% 3) photoactivated region only, during photoactivation
PAregfileglob = 'FRAP_%03d_*_Bleach_*C0.tif';
% 4) general ER luminal marker, during photoactivation
ERfileglob = 'FRAP_%03d_*_Bleach_*C3.tif';
% 5) metadatafile exported from ImageJ (used to get frame rate)
metafile='Metadata_210423_COS7_ER-PAGFP_ER-mCherry_WT.csv';

imgfileglobs = {PAfileglob,prefileglob,PAregfileglob,ERfileglob};

%% Go through, create all cell objects, load in info, manually identify segmentation regions
% Does not do detailed analysis yet
clear allcells
for pcc = 1:length(processnums)
    cc = processnums(pcc); % process this cell / imaging run
    
    % establish specific file names containing img data
    for fc = 1:length(imgfileglobs)
        fname = sprintf([dirname genfilename imgfileglobs{fc}],cc);
        tmp = dir(fname);
        if (isempty(tmp))
            error(sprintf('Failed to find file: %s',fname))
        end
        imgfiles{fc} = tmp.name;
    end
    %%
    % create a cell object and load data into it
    % Note: this does NOT load the actual raw image files (except for one frame
    % saved as CL.ERimg saved for quick visualization
    cellname = sprintf(cellnamepattern,cc);   
    CL = CellObjPA(cellname);
    CL.loadCellData(dirname,imgfiles{1},imgfiles{2},imgfiles{3},imgfiles{4},metafile);
    
    % --------------
    % manually identify rectangular activation region
    CL.getActROI(struct('dodisplay',1,'manualrect',1));
    
    % --------------------
    % manually identify cell and nucleus outlines, based on the general ER region    
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
    
    allcells(pcc) = CL;
end

%% save all cell objects
save(savefile,'allcells')

% ---------
%% Go through and compute wedge ROIs at different distances
% also extract fluorescent signal in each ROI

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

% get non-photoactivated luminal marker traces?
getnonPATrace = false;

%clear regionTraces
for cc = 1:length(allcells)
    %figure(cc)
    CL = allcells(cc);
    
    % Rvals are the radii at the center of each ring band
    [Rvals,whichrad] = CL.getWedgeROIs(maxR,dR,Rwidth,optionsROIs);
    
    CL.cellROI = CL.fullcellROI;
    
    % regionTraces saves traces of fluorescence intensity per area for each
    % wedge or ring ROI
    [regionTraces{cc},~] = CL.getROItraces(getnonPATrace);
    
    %% Plot example ROIs and traces to make sure everything makes sense
    
    % angle (in degrees, btwn -180 and 180) around the circle, indicating roughly
    % which wedge to show for each ring
    % 0 degrees points east. Angle increases clockwise.
    pickangle = -90;
    
    CL.plotExampleWedgeROIs(pickangle)
end

%% save all cell objects, with defined ROIs and signals
save(savefile,'allcells','Rvals','optionsROIs')

%% Go through and fit a double-exponential function to each ROI signal in each cell
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

% ignore all half-times that are more than this factor * total imaging time
optionsFit.maxtscl = 5;

% avghalftime is the median from all cells together
[avghalftime,allhalftimes,allcellind,avghalftimecells] = getAvgHalfTimes(allcells,optionsFit);

%% bootstrap at the cell-to-cell level to get error bars
[medhalftime,stehalftime] = bootstrapwedges(allcells,Rvals,allcellind,allhalftimes);

%% save all cell objects, with defined ROIs and signals and halftimes
save(savefile,'allcells','Rvals','optionsROIs','optionsFit','avghalftime','allhalftimes','allcellind','avghalftimecells')

%% plot results (half times vs Rvals for individual cells and overall median
figure
cmat = jet(length(allcells));
for cc = 1:length(allcells)
    loglog(Rvals,avghalftimecells(:,cc),'.-','Color',cmat(cc,:))
    hold all
end
loglog(Rvals,medhalftime,'k.-','LineWidth',2)
errorbar(Rvals,medhalftime,stehalftime,'k')
hold off
xlim([4,15])
set(gca,'FontSize',14)
xlabel('distance (um)')
ylabel('signal rise time')
title('Many Cells')
