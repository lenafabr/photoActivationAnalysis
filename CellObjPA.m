% class definition for a cell object, with dynamic fluorescent data
% to use for analyzing photoactivation movies

classdef CellObjPA < handle
    
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
            % WARNING: assumes same frame rate for both pre and during
            % photoactivation period
            % -------
            if (exist([dirname metadatafile],'file'))
                data = readtable([dirname metadatafile]);
                
                for kc = 1:size(data,1)
                    if (contains(data.Key{kc},'mage|ATLConfocalSettingDefinition|CycleTime')...
                            & ~contains(data.Key{kc},'Pre Series'))
                        CL.dt = str2num(data.Value{kc});
                        break
                    end
                end
            else
                warning('Failed to load metadata file: %s', [dirname metadatafile])
            end
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
            
            if (exist('options','var'))
                % copy over passed options
                opt = copyStruct(options,opt);
            end
           
            %%
            imPAreg = imread([CL.DirName CL.PAregfile],2);
            
            % Identify centroid and size of photoactivated region
            if (opt.manualrect)
                imshow(imPAreg,[])
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
                %threshold and fill photoactivated region
                T = graythresh(imPAreg);
                imPAbw = imbinarize(imPAreg,T);
                
                % dilate and find connected component to get PA region
                se = strel('disk',opt.dilaterad);
                imPAdil = imdilate(imPAbw,se);
                
                % get largest connected component
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
            CL.actROI.rad = max([PArect(1)+PArect(3)-PAcent(1), PAcent(1)-PArect(1),...
                PAcent(2)-PArect(2), PArect(2)+PArect(4)-PAcent(2)]);
            
            CL.actROI.bound = [PArect(1),PArect(2); PArect(1)+PArect(3),PArect(2); ...
                PArect(1)+PArect(3), PArect(2)+PArect(4); PArect(1), PArect(2)+PArect(4)];
            %%
             if (opt.dodisplay>0)
                imshow(imPAbw)
                hold on
                rectangle('Position',CL.actROI.rect,'Edgecolor','m','LineWidth',2)
                plot(PAcent(1),PAcent(2),'m*')
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
        
        function showCellROI(CL,adjust)            
            % redraw current cell ROI
            % if adjust is present and true, allow for adjusting            
            imshow(CL.ERimg,[])
            
            CL.cellROI.ROI = drawpolygon('Position',CL.cellROI.bound);
            actroi = drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','y');
            if (isempty(CL.actROI.mask))
                CL.actROI.mask = createMask(actroi);
            end
            
            if (exist('adjust','var'))
                if (adjust)
                    input('Hit enter when done adjusting\n')
            
                    CL.cellROI.bound = CL.cellROI.ROI.Position;   
                    CL.cellROI.mask = createMask(CL.cellROI.ROI);
                end
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
            opt.minarea = 0.8 ;        
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
            
            opt.dodisplay = 0;
            
            if (exist('options','var'))
                opt = copyStruct(options,opt);
            end
            
            thetaall = linspace(0,2*pi,opt.ntheta+1);
            thetain = thetaall(1:end-1);
            thetaout = thetaall(2:end);
            
            % mask of cell without nucleus
            nonucmask = CL.fullcellROI.mask & ~CL.nucROI.mask;
            
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
            
            %% get angular slice masks
            if (isnan(opt.arclen))
                slicecoords = zeros(4,2);
                slicemasks = zeros(size(CL.ERimg,1),size(CL.ERimg,2),length(thetain));
                for tc = 1:length(thetain)
                    slicecoords(1,:) = CL.fullcellROI.cent;
                    slicecoords(2,:) = CL.fullcellROI.cent+CL.fullcellROI.rad*[cos(thetain(tc)) sin(thetain(tc))];
                    slicecoords(3,:) = CL.fullcellROI.cent+CL.fullcellROI.rad*[cos(thetaout(tc)) sin(thetaout(tc))];
                    slicecoords(4,:) = CL.fullcellROI.cent;
                    
                    sliceroi = drawpolygon('Position',slicecoords,'Visible','off');
                    slicemasks(:,:,tc) = createMask(sliceroi);
                end
            end
            
            %%
            ct = 0; % count new ROIs
            clear allROIs
            if (opt.dodisplay)
                imshowpair(CL.ERimg,nonucmask)
                cmap = jet(length(Rvalsin));
            end
            
            ringROIs = struct('cent',{},'rad',{},'mask',{},'type',{},'whichrad',{},'whichth',{},'bound',{});
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
                    ringbounds = bwboundaries(ringmask);
                    ringROIs(rc) = struct('cent',CL.actROI.cent,'rad',Rvals(rc)*CL.pxperum,'mask',ringmask2,...
                        'type','ring','whichrad',rc,'bound',NaN,'whichth',NaN);                      
                    ringROIs(rc).bound = ringbounds;
                end
                
                if (~isnan(opt.arclen))
                    % get slices of a fixed arc length
                    totarc = 2*pi*Rvals(rc);
                    
                    if (~isnan(opt.arcshift))
                        % shift wedges by this arc length
                        arcstart = 0:opt.arcshift:totarc-opt.arclen/2;
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
                        slicecoords(1,:) = CL.fullcellROI.cent;
                        slicecoords(2,:) = CL.fullcellROI.cent+CL.fullcellROI.rad*[cos(thetain(tc)) sin(thetain(tc))];
                        slicecoords(3,:) = CL.fullcellROI.cent+CL.fullcellROI.rad*[cos(thetaout(tc)) sin(thetaout(tc))];
                        slicecoords(4,:) = CL.fullcellROI.cent;
                        
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
                    if nnz(wedgemask) > opt.minarea*nw
                        wedgemask = wedgemask & nonucmask;
                        ct = ct+1;
                        wedgebound = bwboundaries(wedgemask);
                        wedgebound = fliplr(wedgebound{1});
                        info = polygeom(wedgebound(:,1),wedgebound(:,2));
                        newROI = struct('mask',wedgemask,'bound',wedgebound,'cent',info(2:3),'rad',Rvals(rc)*CL.pxperum,...
                            'whichrad',rc,'whichth',tc,'type','wedge');
                        whichrad(ct) = rc;
                        allROIs(ct) = newROI;
                        
                        if (opt.dodisplay>1)
                            drawpolygon('Position',allROIs(ct).bound,'Color',cmap(rc,:));
                            drawnow
                        end
                    end
                end
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
        
        function [regionTraces,imgs] = getROItraces(CL,getnonPATrace)
            % get time-traces for all the cell ROIs
            % reads in all the images but does not save them to cell object
            % to conserve space
            % getnonPATrace = also keep track of trace in
            % non-photoactivated signal
            
            if (~exist('getnonPATrace'))
                getnonPATrace = 0;
            end
            
            
            % get a second photoactivated signal as well.
            getsignal2 = ~isempty(CL.PA2file); 
            
            % load in all images
            imgs = loadImages(CL.DirName,CL.PAprefile,CL.PAfile);
            
            if (getsignal2)
                imgs2 = loadImages(CL.DirName,CL.PA2prefile,CL.PA2file);
            end
            if (getnonPATrace)
                imgsER = loadImages(CL.DirName,CL.ERprefile,CL.ERfile);
            end
            
            % put together all region masks
            allmasks = cat(3,CL.ROIs.mask);
            
            % add on activation mask and whole cell mask
            allmasks(:,:,end+1) = CL.actROI.mask;
            allmasks(:,:,end+1) = CL.cellROI.mask;
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
                            
    end
    
end