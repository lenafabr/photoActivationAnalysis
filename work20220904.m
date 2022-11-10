cd ~/UCSD/proj/ERphotoact/photoActivationAnalysis/
%
addpath(genpath('./'))
%% load compiled puff data from Cecile

dirname = '../../ERCaSims/results/capuffs/';
fglob = 'allpuff*_centered.mat';

files = dir([dirname fglob]);
fnames = {files.name};

for fc = 1:length(fnames)
    fname = [dirname fnames{fc}] ;
    load(fname)
    
    tmp = split(fnames{fc},'.mat');
    name = tmp{1}; %variable name
    
    allpuffs{fc} = eval(name);
    namelist{fc} = name;
    
    tmp = regexp(namelist{fc},'allpuff_(\w+)_([0-9]+)_','tokens');
    groupname{fc} = tmp{1}{1};
    dates{fc} = tmp{1}{2};    
    
    if ismember(dates{fc},{'140622','060722'})
        frametime(fc) = 0.66;
    else
        frametime(fc) = 0.15;
    end
end
pufflist = horzcat(allpuffs{:});

%% Load image data
%% general directory name
upperdirname = '/media/ekoslover/PortableDrive2/pc1/.snapshots/snapshot.0/data/proj/ERtransport/calcium_release/';
%gendirname = '/data/proj/ERtransport/calcium_release/COS7-WT-GCaMP3-cIP3-Flash-210422/';
gendirname = [upperdirname 'COS7-WT-GCaMP3-cIP3-Flash-210422/'];

% cell name
cellname = '210422-016-puff';

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


%% test gui
app = roiSignalProcess('CL',CL)

%% save app state and ROIs to file
appstate = app.saveAppState();
CL = app.CL;
appstate.CL.imgs = []; % clear images to avoid large file size
dataname = ['RTN3_', strrep(CL.Name,'-','_'), '_filt'];
assignin('base',dataname,app.CL.ROIs)
%
dirname = '../results/';
fname = [dirname dataname '.mat'];

save(fname,'appstate',dataname)

%% load a previous set of data and restore app
load('../results/WT_210422_010c_puff_filt.mat')
CL.imgs = loadImages(CL.DirName,CL.PAprefile,CL.PAfile);
app.loadAppState(appstate)
app.CL.imgs = CL.imgs;

%% plot results from Cecile
ROIs = allpuff_WT_210422_16c_centered

tvals = 1:length(ROIs(1).avgsignal);
cmat = jet(length(ROIs));
for rc = 1:length(ROIs)
    plot(tvals,ROIs(rc).avgsignal,'Color',cmat(rc,:))
    hold all
end
hold off

%% put Cecile results in app
%app.CL.ROIs = allpuff_WT_210422_3c_centered;

app.CL.ROIs = allpuff_RTN3_010622_10c_centered;

%% ---------- OLD STUFF ----------------------
%% save app state
appstate = struct();

props = properties(app);
values = cell(1, length(props));

for i = 1:length(props)
    propName = props{i}
    property = app.(propName);
    if (strcmp(propName,'img'))
        values{i} = app.img;
    elseif (strcmp(propName,'CL'))
        values{i}= app.CL;
    elseif isprop(property, 'Value')
        values{i} = app.(propName).Value;
    end    
end
appstate.props = props;
appstate.values = values;
appstate.CL.inmgs = [];


%% restore appstate
props = appstate.props;
values = appstate.values;
for i = 1:length(props)
    propName = props{i};
    property = app.(propName);
    if (strcmp(propName,'img'))
        app.img = values{i};
    elseif (strcmp(propName,'CL'))
        app.CL = values{i};
    elseif isprop(property, 'Value')
        app.(propName).Value = values{i};
    end
end

%% save the current state of the app
fields= properties(app);
appstate = struct();
% clear large variable
CL.imgs = [];
for fc = 1:length(fields)
    %fields{fc}
    propName = fields{fc};
    val = app.(propName);
    if (~isa(val,'handle'))
        % save everything that is not a handle;
        appstate.(propName) = val;
    elseif (isprop(val,'Value'))
        % save everything that has a value        
        appstate.(propName) = app.(propName).Value;
    end
end

%% Load app state from saved values
fields = fieldnames(appstate);
for fc = 1:length(fields)
    propName = fields{fc};
    val = appstate.(propName);
    appval = app.(propName);
    if (~isa(appval,'handle'))        
        app.(propName) = val;
    elseif (isprop(appval,'Value'))
        % save everything that has a value
        app.(propName).Value = val;
    end
end


%%
