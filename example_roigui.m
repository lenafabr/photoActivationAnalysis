% general directory name
gendirname = '/data/proj/ERtransport/calcium_release/COS7-WT-GCaMP3-cIP3-Flash-210422/';

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

%%
imshow(CL.imgs(:,:,1),[0,0.3])


%%
addpath('./roisignalGUI/')
%% test gui
app = roiSignalProcess('CL',CL)