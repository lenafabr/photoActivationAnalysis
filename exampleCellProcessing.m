%% This is an example script for loading data involving a particular cell
% and saving to a cell object for future use
% saves one image frame, ROI definitions, and time traces only
% so should not generate very large data files
% This approach avoids redefining cell boundaries, regions, etc, for cells
% that have already been processed


%% create a cell object and load data
CL = CellObjPA('depleted_FRAP004');

dirname = '/data/proj/ERtransport/Edward013020/extracted_data/';
% before photoactivation
prefile = 'C3-190508_COS7_%s_pre.tif';   
% during photoactivation
PAfile = 'C3-190508_COS7_%s_bleach.tif';
% photoactivation region
PAregfile = 'C1-190508_COS7_%s_bleach.tif';
% general ER marker during photoactivation
ERfile = 'C4-190508_COS7_%s_bleach.tif';
% metadata file
metafile = '190508_COS7_%s_metadata.csv';

CL.loadCellData(dirname,PAfile,prefile,PAregfile,ERfile,metafile);

%%
CL.getActROI(struct('dodisplay',1))

%%
CL.setCellROI();

%% show and/or adjust cell ROI
CL.showCellROI(true);

%% get ring regions
ringwidth = 20; % width of ring (in pixels)
dr = 10; % radial separation between rings. Set to ringwidth to make them non-overlapping
nreg = 10; % number of ring regions

CL.getRingROIs(nreg,dr,ringwidth,struct('dodisplay',1))

%% display a particular ring mask to check everything makes sense
imshowpair(CL.ERimg,CL.ROIs(6).mask)
drawrectangle('Position',CL.actROI.rect,'Color','w')

%% get traces over each region
CL.getROItraces()

%% save cell to file 
savefile = sprintf('../celldata/%s_data.mat',CL.Name);
save(savefile,'CL')

%% load from file (START FROM HERE if you have worked with this cell before)
load('../celldata/depleted_FRAP004_data.mat')
%% plot traces vs time
tvals = (1:CL.NFrame)*CL.dt;
regionTraces = cat(1,CL.ROIs.avgsignal);
cmap = jet(length(CL.ROIs))
set(gcf,'defaultAxesColorOrder',cmap)
plot(tvals,regionTraces)
xlabel('time (sec)')
ylabel('avg signal')