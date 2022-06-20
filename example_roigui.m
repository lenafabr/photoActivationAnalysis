% general directory name
upperdirname = '/media/ekoslover/PortableDrive2/pc1/.snapshots/snapshot.0/data/proj/ERtransport/calcium_release/';
%gendirname = '/data/proj/ERtransport/calcium_release/COS7-WT-GCaMP3-cIP3-Flash-210422/';
gendirname = [upperdirname 'COS7-WT-GCaMP3-cIP3-Flash-210422/'];

% cell name
cellname = '210422-015-puff';

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

%% get integrated peaksignals
allpuffinteg = [app.ROIs.puffinteg]

%% plot gaussian peak overlays to see if width is judged correctly
figure(app.PlotFigure)
hold all
xvals0= linspace(-5,5,100);
gauss0 = normpdf(xvals0)*sqrt(2*pi);

for rc = 1:length(app.ROIs)
    puffind = app.ROIs(rc).puffind;
    signal = app.ROIs(rc).avgsignal;
    tlist = 1:length(app.ROIs(rc).avgsignal);
    
    sigmin = min(signal); sigmax = max(signal);
    basesignal = app.ROIs(rc).basesignal;
    
    for cc= 1:length(puffind)
        base = app.ROIs(rc).basesignal(puffind(cc));
        w = app.ROIs(rc).puffwidths(cc);
        gauss = gauss0*(signal(puffind(cc))-base) + base;
        xvals = xvals0*w/2 + tlist(puffind(cc));        
                
        gaussshift  = (gauss-sigmin)/(sigmax-sigmin)+rc;
        
        baseshift = (basesignal-sigmin)/(sigmax-sigmin)+rc;
        plot(tlist,baseshift,'Color',app.roicmap(rc,:))
        plot(xvals,gaussshift,'--','Color',app.roicmap(rc,:))
    end
end
hold off