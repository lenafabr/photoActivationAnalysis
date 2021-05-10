% class definition for a cell object, with dynamic fluorescent data
% to use for analyzing photoactivation movies

classdef CellObjPA < matlab.mixin.Copyable   
    
    properties
        % Cell Name
        Name = 'default';
        
        % image size
        ImgSize = [];
        
        % number of frames
        NFrame = 0;
        
        % directory where image data is stored
        DirName = '';
        
        % filenames:
        PAfile = ''; % photoactivated channel filename
        PAprefile = ''; % photoactivated channel before activation
        PAregfile = ''; % photoactivated region during activation
        % keep track of a second photoactivated signal if needed
        PA2file = ''; % photoactivated channel filename
        PA2prefile = ''; % photoactivated channel before activation
        PA2regfile = ''; % photoactivated region during activation
        ERfile = ''; % constitutive ER marker during activation
        ERprefile = '';
        
        %time interval (seconds per frame)
        dt = 0;
        
        % spatial resolution (px/um)
        pxperum = 0;
        
        % frame at which photoactivation started
        startPA = 0;               
               
        % single snapshot image of ER in this cell
        % uses up memory but not too much since it's just one snapshot
        ERimg = NaN;
        
        % array of ROI structures
        % each ROI structure should have the fields:
        % boundary, mask, radius, center, avgsignal (signal per area over time)
        ROIs = [];
        
        % activated ROI (same fields)
        actROI = [];
        
        % ROI for the part of the cell to be analyzed
        % this has the field ROI for the actual adjustable ROI object
        cellROI = [];
        
        % ROI for entire cell outline
        fullcellROI = [];
        % ROI for nucleus outline
        nucROI = [];
        
        % function used for fitting signals vs time
        fitfunc = NaN;
        
        % save images if desired
        imgs = [];
    end
    
    methods
        function CL = CellObjPA(name)
            % create a cell object with the given name
            CL.Name = name;
            
            % set up empty ROIs structure
            CL.ROIs = struct('bound',{},'mask',{},'rad',{},'cent',{},'avgsignal',{},'avgsignal2',{});
            
            CL.actROI = struct('bound',[],'mask',[],'rad',[],'cent',[],'avgsignal',[],'rect',[],'avgsignal2',[]);
            
            CL.cellROI = struct('bound',[],'mask',[],'rad',[],'cent',[],'avgsignal',[],'ROI',NaN,'avgsignal2',[]);
        end
        
        function loadCellData(CL,dirname,PAfile,PAprefile,PAregfile,ERfile,metafile)
            % load overall data about a cell movie
            % does *NOT* load the images themselves, since those use a lot of
            % memory
            % saves file names, frame time, px size, number images, etc.
            % for each file name, if it contains %s, replace with cell name
            
            % ----------
            %% save the file names for this cell
            % ----------
            CL.DirName = dirname;
            if (contains(PAfile,'%s'))
                CL.PAfile = sprintf(PAfile,CL.Name);
            else
                CL.PAfile = PAfile;
            end
            
            if (contains(PAprefile,'%s'))
                CL.PAprefile = sprintf(PAprefile,CL.Name);
            else
                CL.PAprefile = PAprefile;
            end
            
            if (contains(PAregfile,'%s'))
                CL.PAregfile = sprintf(PAregfile,CL.Name);
            else
                CL.PAregfile = PAregfile;
            end
            
            if (contains(ERfile,'%s'))
                CL.ERfile = sprintf(ERfile,CL.Name);
            else
                CL.ERfile = ERfile;
            end
            
            if (contains(metafile,'%s'))
                metadatafile = sprintf(metafile,CL.Name);
            else
                metadatafile = metafile;
            end
            
            % ------------
            %% load in metadata about images
            % ------------
            % data during photoactivation period
            if (~exist([dirname CL.PAfile],'file'))
                error('Failed to load photoactivated image file: %s', [dirname CL.PAfile])
            end
            infoPA = imfinfo([dirname CL.PAfile]);
            
            % number of frames
            nPAframe = length(infoPA);
            info1 = infoPA(1);
            % number of px in image
            CL.ImgSize = [info1.Height, info1.Width];
            % px per um scale factor
            CL.pxperum = info1.XResolution;
            
            % data pre-photoactivation
            if (~exist([dirname CL.PAprefile],'file'))
                error('Failed to load pre-photoactivated image file: %s', [dirname CL.PAprefile])
            end
            infopre = imfinfo([dirname CL.PAprefile]);
            npreframe = length(infopre);
            
            CL.NFrame = npreframe + nPAframe;
            CL.startPA = npreframe+1;
            
            % -------
            %% load in frame rate from metadata
            % WARNING: only loads frame rate for during the
            % photoactivation period
            % -------
            
            % get FRAP number from cell name
            numstr = regexp(CL.Name,'FRAP([0-9]+)','tokens');
            frapnum = sprintf('FRAP%s',numstr{1}{1});
            frapnum2 = sprintf('FRAP %s',numstr{1}{1});
            
            if (exist([dirname metadatafile],'file'))
                data = readtable([dirname metadatafile],'Delimiter',',');
                
                for kc = 1:size(data,1)
                    if (contains(data.Key{kc},'mage|ATLConfocalSettingDefinition|CycleTime')...
                            && ~contains(data.Key{kc},'Pre Series') && ...
                            (contains(data.Key{kc},frapnum) || contains(data.Key{kc},frapnum2)))
                        CL.dt = str2num(data.Value{kc});
                        break
                    end
                end
            else
                warning('Failed to load metadata file: %s', [dirname metadatafile])
            end
            
            %% load in one frame to be able to quickly visualize the cell
            CL.ERimg = imread([CL.DirName,CL.ERfile],1);
        end
        
        function getActROI(CL,options)
             % automatically identify photoactivated region from the PAreg
            % image file
            
            % ----- set default options ----
            opt = struct();
            
            opt.dodisplay = 0; % display identified region?
            % how much to dilate to identify connected region
            opt.dilaterad = 5;
            % select activated rectangle manually
            opt.manualrect = 0;
            
            % the options below are mostly used to deal with simulated data
            % treated as image
            
            % input image rather than reading from file
            opt.imPAreg = NaN;            
            % apply median filter of some size before binarizing           
            opt.medfilter = NaN;
            % use a percentile cutoff for thresholding
            opt.threshprctile = NaN;
            
            if (exist('options','var'))
                % copy over passed options
                opt = copyStruct(options,opt);
            end
           
            %%
            if (isnan(opt.imPAreg))
                imPAreg = imread([CL.DirName CL.PAregfile],2);
            else
                imPAreg = opt.imPAreg;
            end
            
            % Identify centroid and size of photoactivated region
            if (opt.manualrect)
                imshow(imPAreg,[])
                title(sprintf('%s: identify activated rectangle',CL.Name),'Interpreter','none')
                disp('Draw rectangle')
                rect = drawrectangle;
                input('Hit enter when done adjusting rectangle');
                mask = createMask(rect);
                imPAregmask = imPAreg;
                imPAregmask(~mask(:)) = 0;
                %threshold and fill photoactivated region
                T = graythresh(imPAregmask);
                imPAbw = imbinarize(imPAregmask,T);
            else % identify activated region automatically
                % median filter
                if (~isnan(opt.medfilter))
                    imPAreg = medfilt2(imPAreg,[opt.medfilter,opt.medfilter]);
                end
                
                %threshold and fill photoactivated region
                if (opt.threshprctile==NaN)
                    T = graythresh(imPAreg); % use otsu thresholding
                else
                    % threshold to some percentile                    
                    T = prctile(imPAreg(:),opt.threshprctile);                    
                end
                imPAbw = imbinarize(imPAreg,T);                                
                
                %% dilate and find connected component to get PA region
                se = strel('disk',opt.dilaterad);
                imPAdil = imdilate(imPAbw,se);
                
                %% get largest connected component
                CC = bwconncomp(imPAdil);
                CCsize = cellfun(@(x) length(x), CC.PixelIdxList);
                [a,b] = max(CCsize);
                CClargepx = CC.PixelIdxList{b};
                % and use it to mask out the photoactivated region
                mask = zeros(size(imPAdil));
                mask(CClargepx) = 1;
                imPAbw = imPAbw.*mask;
            end
            %% get bounding box around PA region
            boxinfo = regionprops(double(imPAbw),'BoundingBox');
            PArect = boxinfo.BoundingBox;
            
            % get centroid of PA region
            centinfo = regionprops(double(imPAbw),imPAreg,'WeightedCentroid');
            PAcent = centinfo.WeightedCentroid;                       
            
            CL.actROI.rect = PArect;
            CL.actROI.cent = PAcent;    
                                    
            CL.actROI.rad = min([PArect(1)+PArect(3)-PAcent(1), PAcent(1)-PArect(1),...
                PAcent(2)-PArect(2), PArect(2)+PArect(4)-PAcent(2)]);
            
            CL.actROI.bound = [PArect(1),PArect(2); PArect(1)+PArect(3),PArect(2); ...
                PArect(1)+PArect(3), PArect(2)+PArect(4); PArect(1), PArect(2)+PArect(4)];
            %%
             if (opt.dodisplay>0)
                %imshow(imPAbw)
                imshowpair(CL.ERimg,imPAreg)
                hold on
                rectangle('Position',CL.actROI.rect,'Edgecolor','c','LineWidth',2)
                plot(PAcent(1),PAcent(2),'c*')
                actroi = drawcircle('Center',CL.actROI.cent,'Radius', CL.actROI.rad,'Color','y');
                CL.actROI.mask = createMask(actroi);
                hold off
            end
        end
        
        function setCellROI(CL)
            % manually set the region corresponding to the cell itself
            % must be a continuous polygon, containing the activation
            % center
            % generally avoid going too far around the nucleus
            % hit enter when done!
            
            if (isnan(CL.ERimg))
                 % load ER image
                CL.ERimg = imread([CL.DirName,CL.ERfile],1);
            end      
            
            imshow(CL.ERimg,[]);
            title(CL.Name)
            actroi = drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','y');
            actroirect = drawrectangle(gca,'Position',CL.actROI.rect,'FaceAlpha',0,'Color','m')
            if (isempty(CL.actROI.mask))
                CL.actROI.mask = createMask(actroi);
            end
            
            disp(sprintf(...
                'Click on image to draw polygon outline for overall cell region to analyze. \n Right-click to finish.'))
            CL.cellROI.ROI = drawpolygon(gca);
            
            input('Drag ROI points to adjust. Hit enter when done adjusting\n')
            
            CL.cellROI.bound = CL.cellROI.ROI.Position;     
            CL.cellROI.mask = createMask(CL.cellROI.ROI);
        end        
        
        function setCellNucROI(CL,options)
            % manually set the region corresponding to the full cell and
            % the nucleus (each a continuous polygon containing activation
            % circle)            
            % hit enter for each when done
            
            opt = struct();
            % set cell or nuc as circle with preset center / radius
            opt.cellcent = NaN;
            opt.cellrad = NaN;
            opt.nuccent = NaN;
            opt.nucrad = NaN;
            opt.nucbound = NaN;
            opt.showPA = false
            
            if (exist('options','var'))
                % copy over passed options
                opt = copyStruct(options,opt);
            end            
            
            if (isnan(CL.ERimg) & ~opt.showPA)
                 % load ER image
                CL.ERimg = imread([CL.DirName,CL.ERfile],1);
            end      
            
            if (opt.showPA)
                img = imread([CL.DirName,CL.PAfile],CL.NFrame-CL.startPA-1);
                imshow(img,[]);
            else
                imshow(CL.ERimg,[]);
            end
            title(sprintf('%s: select boundaries', CL.Name),'Interpreter','none')
            actroi = drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','y');
            actroirect = drawrectangle(gca,'Position',CL.actROI.rect,'FaceAlpha',0,'Color','m')
            if (isempty(CL.actROI.mask))
                CL.actROI.mask = createMask(actroi);
            end
            
            % full cell outline
            if (isnan(opt.cellcent) | isnan(opt.cellrad))
                disp(sprintf(...
                    'Click on image to draw polygon outline for full cell. \n Right-click to finish.'))
                ROI = drawpolygon(gca);
                input('Drag ROI points to adjust. Hit enter when done adjusting\n')
            else
                thvals = linspace(0,2*pi,30)';
                bound = opt.cellcent + opt.cellrad*[cos(thvals) sin(thvals)]
                ROI = drawpolygon(gca,'Position',bound);
            end
            CL.fullcellROI = processPolygonROI(ROI);
            
            % nucleus outline
            if (isnan(opt.nuccent) | isnan(opt.nucrad))
                if (isnan(opt.nucbound)) % set nucleus roi manually
                  disp(sprintf(...
                    'Click on image to draw polygon outline for nucleus. \n Right-click to finish.'))
                ROI = drawpolygon(gca,'Color','g');            
                input('Drag ROI points to adjust. Hit enter when done adjusting\n')
                else % set as polygon 
                    ROI = drawpolygon(gca,'Position',opt.nucbound);
                end
            else
                thvals = linspace(0,2*pi,30)';
                bound = opt.nuccent + opt.nucrad*[cos(thvals) sin(thvals)];
                ROI = drawpolygon(gca,'Position',bound);
            end
            CL.nucROI = processPolygonROI(ROI);            
        end        
        
        function showCellROI(CL,adjust)            
            % redraw current ROIs for activated region, cell outline,
            % nucleus
            % if adjust is present and true, allow for adjusting and resave
                                    
            if (~exist('adjust','var'))
                adjust = false;
            end
            
            imshow(CL.ERimg,[])
            actroi =  drawrectangle('Position',CL.actROI.rect,...
                'Color','m','FaceAlpha',0,'InteractionsAllowed','none');
            circroi =  drawcircle('Radius',CL.actROI.rad,'Center',CL.actROI.cent,...
                'Color','m','FaceAlpha',0,'InteractionsAllowed','none');
                
            if (adjust)
                title(sprintf('%s: adjust segmentation',CL.Name),'Interpreter','none')
                cellroi =  drawpolygon('Position',CL.fullcellROI.bound,...
                    'Color',[0,0.7,1],'FaceAlpha',0,'InteractionsAllowed','reshape');
                nucroi =  drawpolygon('Position',CL.nucROI.bound,...
                    'Color','y','FaceAlpha',0,'InteractionsAllowed','reshape');                              
                input('Hit enter when done adjusting\n')
                CL.fullcellROI = processPolygonROI(cellroi);
                CL.nucROI = processPolygonROI(nucroi);
                
            else
                title(sprintf('%s: segmented regions',CL.Name),'Interpreter','none')
                cellroi =  drawpolygon('Position',CL.fullcellROI.bound,...
                    'Color',[0,0.7,1],'FaceAlpha',0,'InteractionsAllowed','none');
                nucroi =  drawpolygon('Position',CL.nucROI.bound,...
                    'Color','y','FaceAlpha',0,'InteractionsAllowed','none');                
            end

        end
        
        function showRings(CL,imgscl)
            % draw the ring regions
            if (~exist('imgscl'))
                imgscl = [];
            end
            imshow(CL.ERimg,imgscl)
            
            drawrectangle('Position',CL.actROI.rect,'Color','m')
            
            hold all
            cmap = jet(length(CL.ROIs));
            for rc = 1:length(CL.ROIs)
                 plot(CL.ROIs(rc).bound(:,1),CL.ROIs(rc).bound(:,2),'-','LineWidth',2,'Color',cmap(rc,:))                 
            end
            hold off
        end
        
        function [Rvals,whichrad] = getWedgeROIs(CL,maxR,dR,ringwidth,options)
            % get wedge-shaped ROIs for the cell
            % store them in CL.ROIs
            % returns central R values for the rings (in um)
            % ROIs.whichrad = index into Rvals associated with each wedge
            % input:
            % maxR = max ring radius in um
            % dR = offset between rings in um
            % ringwidth = width of rings in um            
            
            opt = struct();
            % at least this fraction of the wedge must be in cell area in
            % order to keep it
            opt.minareafrac = 0.8 ;        
            % or have amount total in um^2
            opt.minarea = 0;
            
            opt.ntheta = 30; % default number of angular slices
            % set min radius value in um; if not provided, use the
            % activation radius
            opt.minR = NaN; 
            % optionally, make slices of constant arc length, potentially
            % overlapping
            opt.arclen = NaN;
            % if set, offset all slices by this arc length, regardless of
            % ntheta
            opt.arcshift = NaN;
            % include ring ROIs as well
            opt.getRingROIs = false;
            
            % mask out activated rectangle from regions as well?
            opt.actmask = true;
            
            % erode mask by this many um, to stay away from boundary
            opt.erodemask = 0;
            
            opt.dodisplay = 0;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
                       
            
            thetaall = linspace(0,2*pi,opt.ntheta);
            thetain = thetaall(1:end-1);
            thetaout = thetaall(2:end);
            
            % mask of cell without nucleus
            nonucmask = CL.fullcellROI.mask & ~CL.nucROI.mask;
            % remove activated region as well?
            if (opt.actmask)
                imshow(CL.ERimg,[])
                actrect = drawrectangle('Position',CL.actROI.rect);
                actmask = createMask(actrect);
                nonucmask = nonucmask & ~actmask;
            end
            if (opt.erodemask > 0)
                SE = strel('disk',ceil(opt.erodemask*CL.pxperum),8);
                nonucmask = imerode(nonucmask,SE);
            end
            
            % activation radius in um
            if (isnan(opt.minR))
                actRum = CL.actROI.rad/CL.pxperum;
            else
                actRum = opt.minR;
            end
            
            % inner and outer radii
            Rvalsin = actRum:dR:maxR-dR;
            Rvals = Rvalsin+ringwidth/2;
            Rvalsout = Rvalsin+ringwidth;
            nR = length(Rvalsin);
            
            imshow(CL.ERimg,[])
            title(sprintf('%s: getting wedges',CL.Name),'Interpreter','none')
            
            %% get angular slice masks
            actcent = CL.actROI.cent;
            dists = CL.fullcellROI.bound - actcent;
            dists = sqrt(sum(dists.^2,2));
            cellrad = max(dists);
            if (isnan(opt.arclen))
                slicecoords = zeros(4,2);
                slicemasks = zeros(size(CL.ERimg,1),size(CL.ERimg,2),length(thetain));
                for tc = 1:length(thetain)
                    slicecoords(1,:) = actcent;
                    slicecoords(2,:) = actcent+cellrad*[cos(thetain(tc)) sin(thetain(tc))];
                    slicecoords(3,:) = actcent+cellrad*[cos(thetaout(tc)) sin(thetaout(tc))];
                    slicecoords(4,:) = actcent;
                    
                    sliceroi = drawpolygon('Position',slicecoords,'Visible','on');
                    slicemasks(:,:,tc) = createMask(sliceroi);
                end
            end
            
            %%
            ct = 0; % count new ROIs
            clear allROIs
            if (opt.dodisplay)
                imshowpair(CL.ERimg,nonucmask)
                title(sprintf('%s: getting wedges',CL.Name),'Interpreter','none')
                cmap = jet(length(Rvalsin));
                drawcircle('Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','w')
            end
            
            ringROIs = struct('cent',{},'rad',{},'mask',{},'type',{},'whichrad',{},'whichth',{},'bound',{},'angle',{});
            for rc = 1:nR
                %rc
                outroi = drawcircle('Center',CL.actROI.cent,'Radius',Rvalsout(rc)*CL.pxperum,'Visible','off');
                inroi = drawcircle('Center',CL.actROI.cent,'Radius',Rvalsin(rc)*CL.pxperum,'Visible','off');
                 %outroi = drawcircle('Center',CL.actROI.cent,'Radius',Rvalsout(rc),'Visible','off');
                %inroi = drawcircle('Center',CL.actROI.cent,'Radius',Rvalsin(rc),'Visible','off');
                
                outmask = createMask(outroi);
                inmask = createMask(inroi);
                
                ringmask = outmask & ~inmask;
                
                if (opt.getRingROIs)
                    % save the ring regions as well
                    ringmask2 = ringmask & nonucmask;
                    ringbounds = bwboundaries(ringmask2);
                    for bc = 1:length(ringbounds)
                        ringbounds{bc} = fliplr(ringbounds{bc});
                    end
                    
                    ringROIs(rc) = struct('cent',CL.actROI.cent,'rad',Rvals(rc)*CL.pxperum,'mask',ringmask2,...
                        'type','ring','whichrad',rc,'bound',NaN,'whichth',NaN,'angle',NaN);                      
                    ringROIs(rc).bound = ringbounds;
                end
                
                if (~isnan(opt.arclen))
                    % get slices of a fixed arc length
                    totarc = 2*pi*Rvals(rc);
                    
                    if (~isnan(opt.arcshift))
                        % shift wedges by this arc length
                        arcstart = 0:opt.arcshift:totarc;%-opt.arclen/2;
                    else
                        % fixed total number of wedges
                        arcstart = linspace(0,totarc-opt.arclen/2,opt.ntheta);
                    end
                    
                    arcend = arcstart+opt.arclen;
                    thetain = arcstart/Rvals(rc);
                    thetaout = arcend/Rvals(rc);
                    
                    slicecoords = zeros(4,2);
                    slicemasks = zeros(size(CL.ERimg,1),size(CL.ERimg,2),length(thetain));                   
                    for tc = 1:length(thetain)
                        slicecoords(1,:) = actcent;
                        slicecoords(2,:) = actcent+cellrad*[cos(thetain(tc)) sin(thetain(tc))];
                        slicecoords(3,:) = actcent+cellrad*[cos(thetaout(tc)) sin(thetaout(tc))];
                        slicecoords(4,:) = actcent;
                        
                        sliceroi = drawpolygon('Position',slicecoords,'Visible','off');
                        slicemasks(:,:,tc) = createMask(sliceroi);
                    end
                end
                
                for tc = 1:length(thetain)
                    % intersect with angular slice
                    wedgemask = ringmask & slicemasks(:,:,tc);
                    nw = nnz(wedgemask);
                    
                    % intersect with cell region
                    wedgemask = wedgemask & nonucmask;
                    %imshowpair(CL.ERimg,wedgemask)
                    
                    % only keep if more than 80% within cell region                    
                    
                    nwmask = nnz(wedgemask);
                    
                    if nwmask > opt.minareafrac*nw & nwmask > opt.minarea*CL.pxperum^2
                        wedgemask = wedgemask & nonucmask;
                        ct = ct+1;
                        wedgebound = bwboundaries(wedgemask);
                        boundlen = cellfun(@(x) length(x), wedgebound);
                        [~,bi] = max(boundlen); 
                        wedgebound = fliplr(wedgebound{bi});
                        % suppress polyshape repair warnings
                        warning('off','MATLAB:polyshape:repairedBySimplify')
                        polywedge = polyshape(wedgebound(1:end-1,1),wedgebound(1:end-1,2));
                        [xc,yc] = centroid(polywedge);
                        newROI = struct('mask',wedgemask,'bound',wedgebound,'cent',[xc yc],'rad',Rvals(rc)*CL.pxperum,...
                            'whichrad',rc,'whichth',tc,'type','wedge');
                        
                        % angle associated with wedge center
                        newROI.angle = atan2(yc-actcent(2),xc-actcent(1));
                        
                        whichrad(ct) = rc;
                        allROIs(ct) = newROI;
                        
                        if (opt.dodisplay>1)
                            %drawpolygon('Position',allROIs(ct).bound,'Color',cmap(rc,:));
                            bound = allROIs(ct).bound;
                            hold all
                            plot(bound(:,1),bound(:,2),'LineWidth',2,'Color',cmap(rc,:))
                            drawnow
                        end
                    end
                end
                if (opt.dodisplay>1); hold off; end
            end
            if (opt.getRingROIs)
                % include ring regions
                CL.ROIs = [allROIs ringROIs];
            else
                CL.ROIs = allROIs;
            end
        end
        
        function getRandomROIsRads(CL,cellmaskB,Rvals,regsize,nregperR,options)
            % get random ROIs at specified radiuses around the activation center of the
            % cell
            % Rvals = radius values
            % cellmaskB = cell array of boundaries around masked regions of the cell
            % regsize = size of each region
            % nregperR = number of regions for each R value
            % the ROIs are reasonably well but not perfectly spaced out (subsampling without replacement)
            % and restricted to within the full cell mask minus the nucleus
            
            % sample from this many angular points, as a factor of number of regions
            % desired
            opt.sampleth = 10;
            % number of angular values to use in defining ROI boundaries
            opt.nthregbound = 30;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            allth = linspace(0,2*pi,opt.sampleth*nregperR+1)'; % angles to check
            allth = allth(1:end-1);
            ct = 0;
            clear allROIs whichrad
            
            imshow(CL.ERimg,[])
            
            nR = length(Rvals);
            
            for rc = 1:nR
                % positions along the circle
                Rpts = Rvals(rc)*CL.pxperum*[cos(allth) sin(allth)] + CL.actROI.cent;
                
                % keep only positions within the (possibly disjointed) cell mask
                incellmask = false(size(Rpts,1),1);
                for bc = 1:length(cellmaskB)
                    incellmask = incellmask | inpolygon(Rpts(:,1),Rpts(:,2),cellmaskB{bc}(:,2),cellmaskB{bc}(:,1));
                end
                
                Rpts = Rpts(incellmask,:);
                
                % subsample (without replacement) the desired number of points from
                % within the cell mask at this radius
                indsamp = randsample(nnz(incellmask),nregperR,false);
                Rpts = Rpts(indsamp,:);
                
                for pc = 1:nregperR
                    % get ROIs
                    regroi = drawcircle('Center',Rpts(pc,:),'Radius',regsize*CL.pxperum,'Visible','off');
                    mask = createMask(regroi);
                    
                    newROI = struct('cent',Rpts(pc,:),'rad',regsize*CL.pxperum,'avgsignal',[],'avgsignal2',[],'avgsignalER',[],'whichrad',NaN);
                    newROI.mask = mask;
                    
                    thregbound= linspace(0,2*pi,opt.nthregbound)';
                    newROI.bound = newROI.rad*[cos(thregbound) sin(thregbound)] + newROI.cent;
                    
                    ct = ct+1;
                    allROIs(ct) = newROI;
                    % which radius does each ROI belong to
                    allROIs(ct).whichrad = rc;
                    
                    %disp([rc Rvals(rc)*CL.pxperum ct norm(allROIs(ct).cent-CL.actROI.cent)])
                end
            end
            CL.ROIs = allROIs;
        end

        function getRingROIs(CL,nreg,dr,ringwidth,options)
            % get masks for annular ring regions
            % does NOT set ROI boundaries
            % nreg = number of regions
            % dr = how much radius is shifted for each ring
            % ringwidth = ring width
            
            % ----- set default options ----
            opt = struct();
            
            opt.dodisplay = 0; % display ring outer boundaries?
            opt.nth = 30; % how many angular points to define boundaries
            opt.Rmin = NaN;
            
            if (exist('options','var'))
                % copy over passed options
                opt = copyStruct(options,opt);
            end
                         
            
            if (opt.dodisplay>0)
                if (isnan(CL.ERimg))
                    % load ER image
                    CL.ERimg = imread([CL.DirName,CL.ERfile],1);
                end
                imshow(CL.ERimg,[])
                cmap = jet(nreg);
                drawrectangle(gca,'Position',CL.actROI.rect,'Color','m'); 
                drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','w');                     
                hold all
            end
            
            % first ring starts at approx radius of activated region            
            cent = CL.actROI.cent;
            thlist = linspace(0,2*pi,opt.nth)';        
            if (isnan(opt.Rmin))
                rmin = CL.actROI.rad; % start at activated region radius
            else
                rmin = opt.Rmin; % start at specified minimum radius
            end
            
            for rc = 1:nreg                
                innerrad = rmin+(rc-1)*dr;
                outerrad = rmin + (rc-1)*dr+ringwidth;
                CL.ROIs(rc).rad = (outerrad + innerrad)/2;
                CL.ROIs(rc).cent =cent;
                outerROI = drawcircle(gca,'Center',cent,'Radius',outerrad,'Visible','off');    
                innerROI = drawcircle(gca,'Center',cent,'Radius',innerrad,'Visible','off');
                CL.ROIs(rc).mask = createMask(outerROI);
                innermask = createMask(innerROI);
                CL.ROIs(rc).mask(logical(innermask)) = 0;
                % intersect with cell region
                CL.ROIs(rc).mask = CL.ROIs(rc).mask.*CL.cellROI.mask; 
                
                % set boundary to outercircle
                CL.ROIs(rc).bound = outerrad*[cos(thlist) sin(thlist)] + cent;
                innerbound =  innerrad*[cos(thlist) sin(thlist)] + cent;
                
                
                if (opt.dodisplay>0)                    
                    plot(CL.ROIs(rc).bound(:,1),CL.ROIs(rc).bound(:,2),'-','LineWidth',2,'Color',cmap(rc,:))
                    plot(innerbound(:,1),innerbound(:,2),'--','LineWidth',2,'Color',cmap(rc,:))
                end
            end
            
            if (opt.dodisplay>0)
                hold off
            end                        
        end
        
        function [regionTraces,imgs] = getROItraces(CL,getnonPATrace,loadoptions,options)
            % get time-traces for all the cell ROIs
            % reads in all the images but does not save them to cell object
            % to conserve space
            % getnonPATrace = also keep track of trace in
            % non-photoactivated signal
            
            if (~exist('getnonPATrace'))
                getnonPATrace = 0;                
            end
            if (~exist('loadoptions'))
                loadoptions = struct();
            end
            
            % other options, with defaults
            opt = struct();
            % input images rather than loading from file
            opt.imgs = NaN;
                        
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            if (~exist('loadoptions'))
                loadoptions = struct();
            end
            
            
            % get a second photoactivated signal as well.
            getsignal2 = ~isempty(CL.PA2file); 
            
            % load in all images
            if (isnan(opt.imgs))
                imgs = loadImages(CL.DirName,CL.PAprefile,CL.PAfile,loadoptions);
            else
                imgs = opt.imgs;
            end
            
            if (getsignal2)
                imgs2 = loadImages(CL.DirName,CL.PA2prefile,CL.PA2file,loadoptions);
            end
            if (getnonPATrace)
                imgsER = loadImages(CL.DirName,CL.ERprefile,CL.ERfile,loadoptions);
            end
            
            % put together all region masks
            allmasks = cat(3,CL.ROIs.mask);
            
            % add on activation mask and whole cell mask
            allmasks(:,:,end+1) = CL.actROI.mask;
            % mask for cell region around activation zone; may be undefined
            if (isempty(CL.cellROI.mask)) 
                allmasks(:,:,end+1) = CL.fullcellROI.mask;
            else
                allmasks(:,:,end+1) = CL.cellROI.mask;
            end
            allmasks(:,:,end+1) = CL.fullcellROI.mask;
             
            % area in each mask (in um^2)
            areas = squeeze(sum(sum(allmasks,1),2))/CL.pxperum^2;
            
            % reshape each image into a single vector
            imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
            if (getsignal2)
                img2shape = reshape(imgs2,size(imgs2,1)*size(imgs2,2),size(imgs2,3));
            end
            if (getnonPATrace)
                imgERshape = reshape(imgsER,size(imgsER,1)*size(imgsER,2),size(imgsER,3));
            end
            maskshape = reshape(allmasks,size(imgs,1)*size(imgs,2),size(allmasks,3));
            
            % get matrix of total brightness for each region, for each frame
            % brightness is scaled by the region area
            %regrads = [regROIs.Radius]';
            regionTraces = (maskshape'*imgshape)./areas;
            
            
            nreg = length(CL.ROIs);
            for rc = 1:nreg
                CL.ROIs(rc).avgsignal = regionTraces(rc,:);
            end
            CL.actROI.avgsignal = regionTraces(nreg+1,:);
            CL.cellROI.avgsignal = regionTraces(nreg+2,:);
            CL.fullcellROI.avgsignal = regionTraces(nreg+3,:);
            
            if (getsignal2)
                regionTraces2 = (maskshape'*img2shape)./areas;
                
                for rc = 1:nreg
                    CL.ROIs(rc).avgsignal2 = regionTraces2(rc,:);
                end
                CL.actROI.avgsignal2 = regionTraces2(nreg+1,:);
            end
            
            if (getnonPATrace)
                regionTracesER = (maskshape'*imgERshape)./areas;
                
                for rc = 1:nreg
                    CL.ROIs(rc).avgsignalER = regionTracesER(rc,:);
                end
                CL.actROI.avgsignalER = regionTracesER(nreg+1,:);
                CL.cellROI.avgsignalER = regionTracesER(nreg+2,:);
                CL.fullcellROI.avgsignalER = regionTracesER(nreg+3,:);
            end
            
        end
        
        function showROI(CL,ind)
            % display a specific set of ROIs
            if (isfield(CL.ROIs,'avgsignal')); subplot(1,2,1); end
            imshow(CL.ERimg,[])
            
            hold all  
            cmat = jet(length(ind));       
            
            for rc = 1:length(ind)
                ROI = CL.ROIs(ind(rc));                
                bounds = ROI.bound;
                if (iscell(bounds)) % multiple boundaries around this region
                    for bc = 1:length(bounds)
                        bound = bounds{bc};
                        plot(bound(:,1),bound(:,2),'Color',cmat(rc,:),'LineWidth',2)
                    end
                else
                    plot(bounds(:,1),bounds(:,2),'Color',cmat(rc,:),'LineWidth',2)
                end
            end
            hold off
            title(sprintf('%s',CL.Name),'Interpreter','none')
            
            if (isfield(CL.ROIs,'avgsignal'))
                % plot avg signal if available
                subplot(1,2,2)
                tvals = (1:CL.NFrame)*CL.dt;
                for rc = 1:length(ind)
                    signal = CL.ROIs(ind(rc)).avgsignal;
                    plot(tvals,signal,'Color',cmat(rc,:))
                    hold all
                end
                hold off
                xlabel('time (s)')
                ylabel('signal')
            end
        end
        
        function newROI = adjustROI(CL,ROIstruct)
            % adjust an ROI info structure, defined by a polygonal boundary
            
            imshow(CL.ERimg,[])
            
            roiobj = drawpolygon('Position',ROIstruct.bound);
            
            input('Hit enter when done adjusting region');
            
            newROI = processPolygonROI(roiobj);
        end
        
        function roiinds = ROIsByPoint(CL,pt)
            % find all ROIs containing a given point
            % rounds point to nearest pixel and uses ROI masks
            ptround = round(pt);
            roiinds = [];
            for rc = 1:length(CL.ROIs)
                if CL.ROIs(rc).mask(ptround(2),ptround(1))
                    roiinds(end+1) = rc;
                end
            end
        end
        function closewedge = closestROICent(CL,pt)
            % find the ROI whose center is closest to the desired point
            
            wedgedist = cat(1,CL.ROIs.cent) - pt;
            wedgedist = sqrt(sum(wedgedist.^2,2));
            [~,closewedge] = min(wedgedist);            
        end
        
        function [halftimes,allcfit,fitfunc] = getHalfTimes(CL,options)
            % calculate half-times, from single or double exponential fit for all ROI regions
            % in this cell
            % saves output in CL.ROIs.halftime                        
            
            opt.maxframe = CL.NFrame;
            % type of fitting function:
            % 1exp = single exponential
            % 2exp = double exponential
            opt.fittype = '1exp';
            
            % solve for time when fitting function reaches this frac of
            % maxium
            opt.cutfrac = 0.5;
            
            
            % minimal ratio of end to start signal to try fitting
            opt.minsignalchange = 0;
            
            % which signal to use 
            % 1 = avgsignal
            % 2 = avgsignal2
            opt.whichsignal = 1;
            
            % do the wedge ROIs and the ring ROIs
            opt.dowedge= true;
            opt.doring = true;
            
            % notify on whether the fitting worked
            opt.fitnotify = 'none';
            
            % how far from end to average when considering whether
            % sufficient signal change
            opt.endsignalavg = 20; % in frames
            
            opt.initsignalavg=NaN;
            
            % window span for smoothing activated region to get max val
            opt.actsmooth = 10; 
            % do not use signal if it starts more than partway up
            % between initial and max values of activated region
            % turned off by default
            opt.maxyshift = inf;
            
            % check for too high a fano factor at the end
            opt.maxendfano=-1;
            opt.endframes=20;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            if (isnan(opt.maxframe))
                opt.maxframe = CL.NFrame;
            end            
            
          % disp(opt)
            
            tvals = (1:opt.maxframe)*CL.dt;
            tfit = tvals(CL.startPA:end);
           % if (CL.startPA == 0)
           %     timeshift = 0;
           % else
                timeshift = tvals(CL.startPA);
           % end
            if (strcmp(opt.fittype,'1exp'))
                fitfunc = @(c,t) c(1)*(1-exp(-(t-timeshift)/c(2)))+c(3);
            elseif (strcmp(opt.fittype,'1expnoshift'))               
                fitfunc = @(c,t) c(1)*(1-exp(-(t-timeshift)/c(2)));    
            elseif (strcmp(opt.fittype,'1expbleach'))               
                fitfunc = @(c,t) c(1)*exp(-(t-timeshift)/c(3)).*(1-exp(-(t-timeshift)/c(2)));    
            elseif (strcmp(opt.fittype,'1exptimeshift'))
                fitfunc = @(c,t) c(1)*(1-exp(-(t-c(4))/c(2)))+c(3);
            elseif(strcmp(opt.fittype,'2exp'))
                fitfunc = @(c,t) c(5)*(1 - c(1)*exp(-(t-timeshift)/c(2)) - c(3)*exp(-(t-timeshift)/c(4)));
            elseif (strcmp(opt.fittype,'2expfixlim'))
                % fix initial and final values from data                                
                fitfunc0 = @(c,t) (1 - c(1)*exp(-(t-timeshift)/c(2)) - c(3)*exp(-(t-timeshift)/c(4))); 
                fitfunc = fitfunc0;
            elseif (strcmp(opt.fittype,'2expfixlimbleach'))
                fitfunc0 = @(c,t) exp(-(t-timeshift)/c(5)).*(1 - c(1)*exp(-(t-timeshift)/c(2)) - c(3)*exp(-(t-timeshift)/c(4))); 
                fitfunc = fitfunc0;
            else
                error(sprintf('opt.fittype must be 1exp or 2exp. %s is invalid', opt.fittype))
            end                
            
            if (strcmp(opt.fittype,'1exp'))
                allcfit = inf*ones(length(CL.ROIs),3);
            elseif (strcmp(opt.fittype,'1expnoshift'))
                allcfit = inf*ones(length(CL.ROIs),2);
            elseif (strcmp(opt.fittype,'1exptimeshift'))
                allcfit = inf*ones(length(CL.ROIs),4);
            elseif (strcmp(opt.fittype,'2exp'))
                allcfit = inf*ones(length(CL.ROIs),5);
            elseif (strcmp(opt.fittype,'2expfixlim'))
                allcfit = inf*ones(length(CL.ROIs),4);
            elseif (strcmp(opt.fittype,'2expfixlimbleach'))
                allcfit = inf*ones(length(CL.ROIs),5);
            else
                allcfit = inf*ones(length(CL.ROIs),5);
            end
            
            for rc = 1:length(CL.ROIs)
                ROI = CL.ROIs(rc);

                if (strcmp(ROI.type,'wedge') && ~opt.dowedge)
                    CL.ROIs(rc).halftime= NaN;
                    continue
                elseif (strcmp(ROI.type,'ring') && ~opt.doring)
                    CL.ROIs(Rc).halftime = NaN;
                    continue
                end
                               
                if (opt.whichsignal==1)
                    signal = CL.ROIs(rc).avgsignal(1:opt.maxframe);
                elseif (opt.whichsignal==2)
                    signal = CL.ROIs(rc).avgsignal2(1:opt.maxframe);
                else
                    error("opt.whichsignal must be 1 or 2")
                end
                
                
                % calculate z score of end signal average vs
                % pre-photoactivation average
                premean = mean(signal(1:CL.startPA+1));
                prestd = std(signal(1:CL.startPA+1))/sqrt(CL.startPA);
                endmean = mean(signal(end-opt.endsignalavg:end));
                zscore = (endmean-premean)/prestd;
                
                if (nnz(~isnan(signal(CL.startPA+1:end))) < 10) % not enough datapoints
                     halftimes(rc) = inf;                   
                    CL.ROIs(rc).halftime = inf;
                    continue
                %elseif (mean(signal(end-10:end))/mean(signal(CL.startPA+1:CL.startPA+10))<opt.minsignalchange)
                %elseif (mean(signal(end-opt.endsignalavg:end))/mean(signal(1:CL.startPA+1))<opt.minsignalchange)
                elseif(zscore < opt.minsignalchange)
                    halftimes(rc) = inf;                   
                    CL.ROIs(rc).halftime = inf;
                    continue
                end
                
                if (CL.startPA ==0)
                    yshift = signal(1);
                else
                    % average over initial period 
                    if (isnan(opt.initsignalavg))
                        endinit = CL.startPA;
                    else
                        endinit = opt.initsignalavg;
                    end
                    yshift = mean(signal(1:endinit));                
                end
                if (strcmp(opt.fittype,'2expfixlim') |strcmp(opt.fittype,'2expfixlimbleach') | strcmp(opt.fittype,'1expbleach') )         
                    % max value in activated region
                    %endvals = mean(CL.actROI.avgsignal(end-10:end));
                    %endvals = prctile(CL.actROI.avgsignal,90);
                    %endvals = max(CL.actROI.avgsignal);
                    
                    actsmooth = smooth(CL.actROI.avgsignal,opt.actsmooth);
                    %if (strcmp(opt.fittype,'2expfixlimbleach') | strcmp(opt.fittype,'1expbleach'))
                    % get max of smoothed signal in activated region                    
                    endvals = max(actsmooth);
                    if (CL.startPA<=1)
                        actyshift = actsmooth(1);
                    else                        
                        actyshift = mean(actsmooth(1:CL.startPA-1));
                    end
                    %else
                        % get smoothed signal at the end                        
                    %    endvals = mean(actsmooth(end-10:end));
                    %end
%                     for rc = 1:length(CL.ROIs)
%                         if (strcmp(CL.ROIs(rc).type,'ring') )
%                             val = mean(CL.ROIs(rc).avgsignal(end-10:end));
%                             if (val>endvals); endvals = val; end
%                         end
%                     end
                end
                
                % check for too much error at the end
                if (opt.maxendfano>0)
                    meanend = mean(signal(end-opt.endframes:end));
                    stdend = std(signal(end-opt.endframes:end));
                    
                    if (stdend/meanend>opt.maxendfano)
                        % too noisy
                        % don't bother trying to fit.
                        CL.ROIs(rc).halftime= NaN;
                        halftimes(rc) = NaN;
                        CL.ROIs(rc).cfit = [];
                        continue
                    end
                end
                
                if (strcmp(opt.fittype,'1exp'))
                    cguess = [max(signal),10,yshift];
                    lb = [0 0 -inf]; ub = [inf inf inf];
                elseif (strcmp(opt.fittype,'1expnoshift'))
                    cguess = [max(signal),10];
                    lb = [0 0]; ub = [inf inf];                   
                    fitfunc = @(c,t) c(1)*(1-exp(-(t-timeshift)/c(2))) + yshift;  
                elseif (strcmp(opt.fittype,'1expbleach'))
                    %cguess = [max(signal),10,100];
                    %lb = [0 0 0]; ub = [inf inf inf];                   
                    %fitfunc = @(c,t) c(1)*exp(-(t-timeshift)/c(3)).*(1-exp(-(t-timeshift)/c(2))) + yshift; 
                    
                    cguess = [endvals,10, 100];
                    lb = [0 0 0]; ub = [inf, inf, inf];
                    fitfunc = @(c,t) c(1)*exp(-(t-timeshift)/c(3)).*(1-exp(-(t-timeshift)/c(2))) + yshift;
    
                elseif (strcmp(opt.fittype,'1exptimeshift'))
                    cguess = [max(signal),10,yshift,0];
                    lb = [0 0 0 0]; ub = [inf inf inf tfit(end)/2];
                elseif (strcmp(opt.fittype,'2expfixlim'))
                     cguess = [0.5 10 0.5 50];       
                     lb = [0 0 0 0]; ub = [1 inf 1 inf];
                     fitfunc = @(c,t) (endvals - actyshift-yshift)*fitfunc0(c,t) + yshift;   
                elseif (strcmp(opt.fittype,'2expfixlimbleach'))
                    
                    cguess = [0.5 1 2 30];
                    lb = [0 0 0 1]; ub = [1 inf 1 inf];
                    maxval = max(CL.actROI.avgsignal);
                    endvals = maxval;
                    fitfunc = @(c,t) maxval*exp(-(t-timeshift)/c(4)).*...
                    (1 - c(1)*exp(-(t-timeshift)/c(2)) - (1-c(1))*exp(-(t-timeshift)/c(3))) + yshift; 

                    % cguess = [0.5 10 0.5 50,10];       
                    % lb = [0 0 0 0 1]; ub = [1 inf 1 inf inf];
                    % fitfunc = @(c,t) (endvals - yshift)*fitfunc0(c,t) + yshift;   
                else
                    cguess = [0.5 10 0.5 100,max(signal)];
                    lb = [0 0 0 0 0]; ub = [1 inf 1 inf inf];    
                end
                
                options=optimoptions('lsqcurvefit','Display',opt.fitnotify,...
                    'MaxFunctionEvaluations',5000, 'MaxIterations',5000,'FunctionTolerance',1e-4);            
                [cfit, resnorm,residual,exitflag,output] = lsqcurvefit(fitfunc,cguess,tfit,signal(CL.startPA:end),lb,ub,options);
                                
                if (exitflag <= 0)
                    CL.ROIs(rc).halftime= NaN;
                    halftimes(rc) = NaN;
                    continue
                end
                
                allcfit(rc,1:length(cfit)) = cfit;                   
                
                if (strcmp(opt.fittype,'1exp') || strcmp(opt.fittype,'1expnoshift') || strcmp(opt.fittype,'1exptimeshift'))                    
                    halftimes(rc) = -log(1-opt.cutfrac)*cfit(2);  
                elseif (strcmp(opt.fittype,'2expfixlim') )     
                    
                    % do not use this wedge if signal starts more than
                    % partway up the activated region signal scale
                    if (yshift > actyshift + opt.maxyshift*(endvals-actyshift))
                         CL.ROIs(rc).halftime= NaN;
                         halftimes(rc) = NaN;
                         continue
                    end
                    
                    cutsignal = opt.cutfrac*(endvals-yshift)+yshift;
                    
                    fvals = fitfunc(cfit,tfit)-cutsignal;
                    nosol = false;                    
                    posind = find(fvals>0);
                    
                    % get left and right boundary for the zero
                    if (isempty(posind))
                        %%
                        tleft = tfit(end);
                        
                        tend = tfit(end);
                        for tc = 1:5
                            tend = tend*2;
                            fend = fitfunc(cfit,tend)-cutsignal;
                            if (fend>0)
                                tright = tend;
                                break
                            end
                        end
                        if (fend<0) % failed to find solution
                            halftimes(rc) = NaN;   
                            nosol= true;
                        end
                    else
                         if posind(1)==1                    
                            warning('function always above cutsignal')
                            CL.ROIs(rc).halftime= NaN;
                            halftimes(rc) = NaN;
                            continue
                         else
                            tleft = tfit(posind(1)-1);
                         end
                        
                        tright = tfit(posind(1));
                    end
                   
                    if (~nosol)
                        [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc(cfit,t)-cutsignal, [tleft tright]);                    
                    end                                                            
                    
%                     [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc(cfit,t)-(opt.cutfrac*(endvals-yshift)+yshift), mean(tfit));                 
%                     if (exitflag<0)
%                         [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc(cfit,t)-(opt.cutfrac*(endvals-yshift)+yshift), max(tfit)*2);                 
%                         if (exitflag<0)
%                             halftimes(rc) = NaN;
%                         end
%                     end
                elseif (strcmp(opt.fittype,'2expfixlimbleach')) 
                    %dfitfunc = @(t) -maxval/cfit(4)*exp(-(t-timeshift)/cfit(4)).*...
                    %(1-cfit(1)*exp(-(t-timeshift)/cfit(2)) - (1-cfit(1))*exp(-(t-timeshift)/cfit(3))) ...
                   % + maxval*exp(-(t-timeshift)/cfit(4)).*...
                    %(cfit(1)/cfit(2)*exp(-(t-timeshift)/cfit(2)) + (1-cfit(1))/cfit(3)*exp(-(t-timeshift)/cfit(3))); 

                   fitfunc2 = @(c,t) maxval.*exp(-(t-timeshift)/c(4)).*...
                   (1 - c(1)*exp(-(t-timeshift)/c(2)) - (1-c(1))*exp(-(t-timeshift)/c(3))) + yshift; 
                
                    %fitfunc2 = @(c,t) maxval.*...
                    %(1 - c(1)*exp(-(t-timeshift)/c(2)) - (1-c(1))*exp(-(t-timeshift)/c(3))) + yshift; 
                
                    % time of max
                    %tmax = fzero(dfitfunc,mean(tfit));
                    %vmax = fitfunc(cfit,tmax);
                    tcheck = linspace(0,max(cfit(2:4))*10,1e4);
                    vmax = max(fitfunc2(cfit,tcheck));
                    %cutsignal = opt.cutfrac*(maxval - yshift)+yshift;
                    cutsignal = opt.cutfrac*(vmax - yshift)+yshift;
                    
                    fvals = fitfunc2(cfit,tfit)-cutsignal;
                    nosol = false;
                    posind = find(fvals>0);
                    
                    % get left and right boundary for the zero
                    if (isempty(posind))
                        %%
                        tleft = tfit(end);
                        
                        tend = tfit(end);
                        for tc = 1:5
                            tend = tend*2;
                            fend = fitfunc2(cfit,tend)-cutsignal;
                            if (fend>0)
                                tright = tend;
                                break
                            end
                        end
                        if (fend<0) % failed to find solution
                            halftimes(rc) = NaN;
                            nosol= true;
                        end
                    else
                        if posind(1)==1
                            warning('function always above cutsignal')
                            CL.ROIs(rc).halftime= NaN;
                            halftimes(rc) = NaN;
                            continue
                        else
                            tleft = tfit(posind(1)-1);
                        end
                        
                        tright = tfit(posind(1));
                    end
                    
                    if (~nosol)
                        [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc2(cfit,t)-cutsignal, [tleft tright]);
                    end
                    roi.fitfunc2 = fitfunc2;
                elseif (strcmp(opt.fittype,'1expbleach'))
                    dfitfunc = @(t) -cfit(1)/cfit(3)*exp(-(t-timeshift)/cfit(3)).*(1-exp(-(t-timeshift)/cfit(2)))...
                        + cfit(1)*exp(-(t-timeshift)/cfit(3)).*1/cfit(2)*exp(-(t-timeshift)/cfit(2));
                    % time of max
                    tcheck = linspace(1e-4,max(cfit(2:3))*10,1e4);
                    vmax = max(fitfunc(cfit,tcheck));
                    
                    %trytimes = linspace(tvals(10),max(cfit(2:3)),50);
                    halfopt = optimset('display','off');
                    % cutsignal = opt.cutfrac*(max(CL.actROI.avgsignal) - yshift)+yshift;
                    cutsignal = opt.cutfrac*(vmax - yshift)+yshift;
                    %for tc = 1:length(trytimes)
                    [halftimes(rc),fval,exitflaghalf] = fzero(@(t) fitfunc(cfit,t)-cutsignal, tfit(2),halfopt);
                    if (exitflaghalf<0) % try different starting points
                        [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc(cfit,t)-cutsignal, tfit(end)*2,halfopt);
                    end
                    
                    %halftimes(rc) = -log(1-opt.cutfrac)*cfit(2)+timeshift;
                else
                    cutval = cfit(5)*(1-cfit(1)-cfit(3)) + opt.cutfrac*(cfit(5)*(cfit(1)+cfit(3)));
                    [halftimes(rc),fval,exitflag] = fzero(@(t) fitfunc(cfit,t)-cutval, mean(tfit));
                end
                
                if (strcmp(opt.fittype,'2expfixlim') | strcmp(opt.fittype,'2expfixlimbleach') | strcmp(opt.fittype,'1expbleach'))
                    halftimes(rc) = halftimes(rc)-timeshift;
                end
                CL.ROIs(rc).halftime = halftimes(rc);
                CL.ROIs(rc).cfit = cfit;
                CL.ROIs(rc).fitfunc = fitfunc;
            end
            %disp('foo')
        end
        
        function plotExampleWedgeROIs(CL,pickangle)
            % plot example wedge ROIs at different radii
            % all lying approximately along the desired angle
            % pickangle = degrees, between -180 and 180 around the
            % photoactivation center
            % 0 degrees corresponds to East, angle increases clockwise
                        
            whichrad = [CL.ROIs.whichrad]; % which radius each roi corresponds to
            % is this a wedge ROI?
            iswedge= cellfun(@(x) contains(x,'wedge'),{CL.ROIs.type});
            % angle for each ROI
            angs = [CL.ROIs.angle];
            
            subplot(1,2,1)
            imshow(CL.ERimg,[]);
            title(sprintf('%s: example ROIs',CL.Name),'Interpreter','none');
            drawrectangle('Position',CL.actROI.rect,'Color','m','FaceAlpha',0.1,'InteractionsAllowed','none')
            drawcircle('Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','m','InteractionsAllowed','none','FaceAlpha',0)
            
            subplot(1,2,2)
            tvals = (1:CL.NFrame)*CL.dt;
            plot(tvals,CL.actROI.avgsignal,'m','LineWidth',1)
            hold all
            
            cmat = parula(max(whichrad));
            for rc  = 1:max(whichrad)
                wedgeind = find(iswedge & whichrad==rc);
                if (isempty(wedgeind)); continue; end
                [~,wcc] = min(abs(angs(wedgeind)-pickangle*pi/180));
                wc = wedgeind(wcc);
                roi = CL.ROIs(wc);
                
                subplot(1,2,1)
                drawpolygon('Position',roi.bound,'Color',cmat(rc,:),'InteractionsAllowed','none','FaceAlpha',0);
                
                subplot(1,2,2)                
                plot(tvals,roi.avgsignal,'Color',cmat(rc,:),'LineWidth',1)
                hold all
                
                % plot fitted function
                if (isfield(roi,'cfit'))
                    if (~isempty(roi.cfit))
                    tfit = tvals(CL.startPA:end);
                    fitvals = roi.fitfunc(roi.cfit,tfit);
                    plot(tfit,fitvals,'--','Color',cmat(rc,:),'LineWidth',1.5)
                    end
                end
            end
            subplot(1,2,2)
            set(gca,'FontSize',14)
            xlabel('time (sec)')
            ylabel('signal per area')
            hold off
        end
        
    end
   
end