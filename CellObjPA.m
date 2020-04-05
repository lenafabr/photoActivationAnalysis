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
        ERfile = ''; % constitutive ER marker during activation
        
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
        
    end
    
    methods
        function CL = CellObjPA(name)
            % create a cell object with the given name
            CL.Name = name;
            
            % set up empty ROIs structure
            CL.ROIs = struct('bound',{},'mask',{},'rad',{},'cent',{},'avgsignal',{});
            
            CL.actROI = struct('bound',[],'mask',[],'rad',[],'cent',[],'avgsignal',[],'rect',[]);
            
            CL.cellROI = struct('bound',[],'mask',[],'rad',[],'cent',[],'avgsignal',[],'ROI',NaN);
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
                error(sprintf('Failed to load photoactivated image file: %s', [dirname CL.PAfile]))
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
                error(sprintf('Failed to load pre-photoactivated image file: %s', [dirname CL.PAprefile]))
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
                    if contains(data.Key{kc},'mage|ATLConfocalSettingDefinition|CycleTime')
                        CL.dt = str2num(data.Value{kc})
                        break
                    end
                end
            else
                warn(sprintf('Failed to load metadata file: %s', [dirname metadatafile]))
            end
        end
        
        function getActROI(CL,options)
             % automatically identify photoactivated region from the PAreg
            % image file
            
            % ----- set default options ----
            opt = struct();
            
            opt.dodisplay = 0; % display identified region?
            
            if (exist('options','var'))
                % copy over passed options
                opt = copyStruct(options,opt);
            end
           
            imPAreg = imread([CL.DirName CL.PAregfile],2);
            
            % Identify centroid and size of photoactivated region
            
            %threshold and fill photoactivated region
            T = graythresh(imPAreg);
            imPAbw = imbinarize(imPAreg,T);
            
            % get bounding box around PA region
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
            
             if (opt.dodisplay>0)
                imshow(imPAbw)
                hold on
                rectangle('Position',CL.actROI.rect,'Edgecolor','m','LineWidth',2)
                plot(PAcent(1),PAcent(2),'m*')
                actroi = drawcircle('Center',CL.actROI.cent,'Radius', CL.actROI.rad,'Color','y')
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
            actROI = drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','y');
            if (isempty(CL.actROI.mask))
                CL.actROI.mask = createMask(CL.actROI);
            end
            
            CL.cellROI.ROI = drawpolygon(gca);
            
            input('Hit enter when done adjusting\n')
            
            CL.cellROI.bound = CL.cellROI.ROI.Position;     
            CL.cellROI.mask = createMask(CL.cellROI.ROI);
        end        
        
        function showCellROI(CL,adjust)            
            % redraw current cell ROI
            % if adjust is present and true, allow for adjusting            
            imshow(CL.ERimg,[])
            
            CL.cellROI.ROI = drawpolygon('Position',CL.cellROI.bound)
            actROI = drawcircle(gca,'Center',CL.actROI.cent,'Radius',CL.actROI.rad,'Color','y');
            if (isempty(CL.actROI.mask))
                CL.actROI.mask = createMask(CL.actROI);
            end
            
            if (exist('adjust','var'))
                if (adjust)
                    input('Hit enter when done adjusting\n')
            
                    CL.cellROI.bound = CL.cellROI.ROI.Position;   
                    CL.cellROI.mask = createMask(CL.cellROI.ROI);
                end
            end
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
            for rc = 1:nreg                
                innerrad = CL.actROI.rad+(rc-1)*dr;
                outerrad = CL.actROI.rad + (rc-1)*dr+ringwidth;
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
        
        function [regionTraces,imgs] = getROItraces(CL)
            % get time-traces for all the cell ROIs
            % reads in all the images but does not save them to cell object
            % to conserve space
            
            % load in all images
            imgs = loadImages(CL.DirName,CL.PAprefile,CL.PAfile);
            
            % put together all region masks
            allmasks = cat(3,CL.ROIs.mask);
            
            % add on activation mask and whole cell mask
            allmasks(:,:,end+1) = CL.actROI.mask;
            allmasks(:,:,end+1) = CL.cellROI.mask;
            
            % area in each mask (in px^2)
            areas = squeeze(sum(sum(allmasks,1),2));
            
            % reshape each image into a single vector
            imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
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
        end
    end
end