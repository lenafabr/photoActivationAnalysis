% for several cell movies, load in information and store in CellObjPA
% object for use in later analysis
% NOTE: this assumes the movies have been saved as .tif files with
% only one channel in each file. 
% It also assumes the pre-bleach and bleach series are stored in separate
% files
% WARNING: for the moment, it uses the frame-time from the bleach series
% for the pre-bleach as well. This is wrong (prebleach frame times are
% different), but the duration of the prebleach period does not currently
% matter for anything.
% Metadata is assumed saved in a .csv file (generated by Fiji when loading
% series from a .lif)

% directory name where all the raw image files are kept
dirname = '/data/proj/ERtransport/Edward021920/extracted_data/';

%% reload previously defined cell boundaries

% if this is set to true then assume you already have a cell object
% and do not need to do the manual segmentation again
% by default, leave as false.
redo = false;
clear allcells
% go through many cell movies
% for multiple cells within a single movie, will need to run it through
% more than once
%%
for cc = 1:2
    
    %% ------------------------------
    % File names where data is stored
    % this is the general first part of the file name
    genfname = '200219_COS7_ER_PA_GFP_RTN4a_';
    % photoactivated channel, before photoactivation
     prefile = [genfname sprintf('FRAP%03d_pre_C2.tif',cc)];
    % photoactivated channel, during photoactivation
     PAfile = [genfname sprintf('FRAP%03d_bleach_C2.tif',cc)];
    % photoactivation region
    PAregfile = [genfname sprintf('FRAP%03d_bleach_C0.tif',cc)];
    % general ER marker during photoactivation
    ERfile = [genfname sprintf('FRAP%03d_bleach_C3.tif',cc)];
    % general ER marker before photoactivation
    ERprefile = [genfname sprintf('FRAP%03d_pre_C3.tif',cc)];      
    % metadata file
    metafile = [genfname sprintf('metadata.csv',cc)];    
    % ---------------------------
    
   
    % give this cell a name so you can recognize it later
    % NOTE: Cell name MUST contain FRAP001 (or some other number)
    name = sprintf('Rtn4a_overexp200219_FRAP%03d',cc);
    
    % data on this cell will be saved in this file
    savefile = sprintf('../celldata/%s_data.mat',name);
    if (redo)
        % reload cell outline from file
        load(savefile,'CL')
        savedcellROI = CL.cellROI;
    end
    
    % create a cell object and load data    
    CL = CellObjPA(name);
    
    CL.loadCellData(dirname,PAfile,prefile,PAregfile,ERfile,metafile);    
    CL.ERprefile = ERprefile;
    
    %% automatically identify activation region
    % if there are multiple rectangles, this will pick only the biggest one
    % if you want to select a region manually (or use other rectangles)
    % then run with ...'manualrect',1).
    CL.getActROI(struct('dodisplay',1,'manualrect',0));    
    %% manually segment the outer boundary of the cell and the nucleus
    if (redo)
        CL.cellROI = savedcellROI;
    else
        CL.setCellNucROI();
    end        
    
    %% get ring and wedge regions
    % Note: rings and wedges can be overlapping 
    % (if dR is smaller than Rwidth, or arcshift is smaller than arclen
    dR = 2; % ring separations (in um)
    Rwidth = 2; % width of rings (in um);
    maxR = 16; % max outer radius (in um)
    minR = 2; % min inner radius (in um)
    ntheta = 30; % number of angular wedges
    arclen = 2; % arc length in um for each wedge
    arcshift = 2; % arc length shift in um between neighbor wedges
    
    
    options= struct('dodisplay',2,'arclen',arclen,'arcshift',arcshift,'minR',minR,'getRingROIs',true);
    [Rvals,whichrad] = CL.getWedgeROIs(maxR,dR,Rwidth,options);
    
    %% get traces over each region 
    % this involves reading in all the image frames and can be a bit slow
    % traces are stored in CL.ROIs.avgtrace
    CL.cellROI = CL.fullcellROI;
    [regionTraces,imgs] = CL.getROItraces(true);
    
    %% save cell to file
    savefile = sprintf('../celldata/%s_data.mat',CL.Name);
    display(sprintf('Saving to file %s', savefile));
    save(savefile,'CL','Rvals')
    
   allcells(cc) = CL;
end

% ---------------
%% Save info for all the cells together in one .mat file
save('../celldata/allcells_Rtn4a_overexp200219.mat','allcells','dR','Rwidth','Rvals','arclen','arcshift','minR','ntheta','maxR')

% -----------------------
%% START FROM HERE if you have loaded the cell before and just want to look at traces
% ------------------------
% load in a previously saved cell data file
load('../celldata/allcells_ctrl200219.mat')
CL = allcells(1);

%% plot traces vs time for each ring
tvals = (1:CL.NFrame)*CL.dt;
regtypes = {CL.ROIs.type};
ringreg = find(strcmp(regtypes,'ring'));
%
cmap = jet(length(ringreg));

subplot(1,2,1)
% show wedge regions
imshow(CL.ERimg,[])
drawcircle('Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','m')
hold all
for rc = 1:length(ringreg)
    bounds = CL.ROIs(ringreg(rc)).bound;
    for bc = 1:length(bounds)
        plot(bounds{bc}(:,1),bounds{bc}(:,2),'LineWidth',2,'Color',cmap(rc,:))
    end
end
hold off
title(sprintf('%s',CL.Name),'Interpreter','none')

subplot(1,2,2)
regionTraces = cat(1,CL.ROIs(ringreg).avgsignal);
%regionTracesER = cat(1,CL.ROIs.avgsignalER);
set(gcf,'defaultAxesColorOrder',cmap)
plot(tvals,regionTraces)
hold all
plot(tvals,CL.actROI.avgsignal,'k')
hold off
xlabel('time (sec)')
ylabel('avg signal')

%% plot traces vs time for each wedge along a particular ring

Rind = 4; % index of the ring to use

regtypes = {CL.ROIs.type};
whichrad = [CL.ROIs.whichrad];
wedgereg = find(strcmp(regtypes,'wedge') & whichrad == Rind);
cmap = jet(length(wedgereg))

subplot(1,2,1)
% show wedge regions
imshow(CL.ERimg,[])
drawcircle('Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','m')
hold all
for wc = 1:length(wedgereg)
    bound = CL.ROIs(wedgereg(wc)).bound    
    plot(bound(:,1),bound(:,2),'LineWidth',2,'Color',cmap(wc,:))
end
hold off
title(sprintf('%s',CL.Name),'Interpreter','none')

subplot(1,2,2)
Rind = 4; % index of the ring to use

regionTraces = cat(1,CL.ROIs(wedgereg).avgsignal);
set(gcf,'defaultAxesColorOrder',cmap)
plot(tvals,regionTraces)
xlabel('time (sec)')
ylabel('avg signal')

% -----------
%% combine all ctrl and all RTN4OE cells (from multiple dates)
% -------------
load('../celldata/allcells_ctrl200213.mat')
allcellsctrl = allcells;
load('../celldata/allcells_ctrl200219.mat')
allcellsctrl = [allcellsctrl allcells];
clear allcells
%
save('../celldata/allcells_ctrl.mat')

%%
load('../celldata/allcells_Rtn4a_overexp200213.mat')
allcellsRTN4OE = allcells;
load('../celldata/allcells_Rtn4a_overexp200219.mat')
allcellsRTN4OE= [allcellsRTN4OE allcells];
clear allcells

save('../celldata/allcells_RTN4OE.mat')