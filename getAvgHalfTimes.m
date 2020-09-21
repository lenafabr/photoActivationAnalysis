function [avghalftime,allhalftimes] = getAvgHalfTimes(allcells,options)
% for a list of cells, calculate average half-time among all regions at a
% given radius
% Rvals lists the radii

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

if (exist('options','var'))
    % copy over passed options
    opt = copyStruct(options,opt);
end

opt

allhalftimes = {};


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
        end
    end
    
    [halftimes,allcfit,fitfunc] = getHalfTimes(CL,opt);
    
    %%
    cmat = jet(nring);
    tfit = tvals(CL.startPA+1:end);
    for rc = 1:nring
        % wedges for this ring index
        % keep only those where half-time is < 2 * movie length
        wedgeind = find(~isring & whichring == rc & halftimes < max(tvals)*opt.maxtscl);
        
        allhalftimes{rc} = [allhalftimes{rc} halftimes(wedgeind)];
        
        if (opt.dodisplay & ~isempty(wedgeind))
            % plot fit for 1 wedge from each radius
            plot(tvals,CL.ROIs(wedgeind(1)).avgsignal,'Color',cmat(rc,:))
            hold all
            plot(tfit,fitfunc(allcfit(wedgeind(1),:),tfit),'--','Color',cmat(rc,:))
        end
    end
    
    if (opt.dodisplay)
        hold off
    end
    
    for rc = 1:length(allhalftimes)
        avghalftime(rc) = mean(allhalftimes{rc});
    end
end
end