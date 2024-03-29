function [avghalftime,allhalftimes,allcellind,avghalftimecells,allwedgeind,plotwedges] = getAvgHalfTimes(allcells,options)
% for a list of cells, calculate average half-time among all regions at a
% given radius
% Rvals lists the radii
% allhalftimes = cell list of half times for each radius value (from all
% cells together)
% allcellind = the index of the cell all the half times in allhalftimes
% correspond to

opt = struct();

% display exponential fits for some regions
opt.dodisplay = 0;

% cut off halftimes above this factor * movie length
opt.maxtscl = 5;

% maximum frame
opt.maxframe = NaN;

% type of fit
opt.fittype = '1exp';

% minimal ratio of end to start signal to try fitting
opt.minsignalchange = 0;

opt.cutfrac = 0.5;

% which signal to process (avgsignal or avgsignal2)
opt.whichsignal=1;


if (exist('options','var'))
    % copy over passed options
    opt = copyStruct(options,opt,'addnew',true);
end

opt

allhalftimes = {};
allcellind = {};
allwedgeind = {};
%%
for cc = 1:length(allcells)
    CL = allcells(cc);
    disp(CL.Name)
    tvals = (1:CL.NFrame)*CL.dt;
    
    isring = strcmp({CL.ROIs.type},'ring');
    whichring = [CL.ROIs.whichrad];
    nring = max(whichring);
    if (nring>length(allhalftimes)) % extend all halftimes to include extra rings
        for rc = (length(allhalftimes)+1):nring
            allhalftimes{rc} = [];
            allcellind{rc} = [];
            allwedgeind{rc} = [];
        end
    end
    
    [halftimes,allcfit,fitfunc] = getHalfTimes(CL,opt);
    
    %%
    cmat = jet(nring);
    tfit = tvals(CL.startPA+1:end);
    plotwedges = [];
    for rc = 1:nring
        % wedges for this ring index
        % keep only those where half-time is < maxtscl * movie length
        wedgeind = find(~isring & whichring == rc & halftimes < max(tvals)*opt.maxtscl);
        
        allhalftimes{rc} = [allhalftimes{rc} halftimes(wedgeind)];
        allcellind{rc} = [allcellind{rc} cc*ones(1,length(wedgeind))];
        allwedgeind{rc} =  [allwedgeind{rc} wedgeind];
        
        if (opt.dodisplay & ~isempty(wedgeind))
            % plot fit for 1 wedge from each radius
            if (opt.whichsignal==1)
                plot(tvals,CL.ROIs(wedgeind(1)).avgsignal,'Color',cmat(rc,:))
            else
                plot(tvals,CL.ROIs(wedgeind(1)).avgsignal2,'Color',cmat(rc,:))
            end
            hold all
            plot(tfit,fitfunc(allcfit(wedgeind(1),:),tfit),'--','Color',cmat(rc,:))

            plotwedges(rc) = wedgeind(1);
        end
        
        % get avg halftimes for this cell only
        avghalftimecells(rc,cc) = nanmedian(halftimes(wedgeind));
    end
    
    if (opt.dodisplay)
        hold off
    end
    
    for rc = 1:length(allhalftimes)
        avghalftime(rc) = nanmedian(allhalftimes{rc});
    end
end
end