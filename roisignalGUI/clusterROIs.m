function [peakregcent,peakregmask] = clusterROIs(ROIs,options)
%% combine together ROIs with similar peak positions and adjacent in space

% returns centers and masks for the clustered regions

opt = struct();
% where meaningful part of signal starts
opt.startPA = 1;

% time cutoff for matching puffs, in terms of frames
opt.timesepcut = 5;

if (exist('options','var'))
    % copy over passed options
    opt = copyStruct(options,opt);
end

%% define a similarity between each pair of regions with puffs
% metric = smallest separation between puff locations
clear signaldist;
tlist = 1:length(ROIs(1).avgsignal(opt.startPA:end));

for rc = 1:length(ROIs)    
    puff = ROIs(rc).puffind;
    signal = ROIs(rc).avgsignal(opt.startPA:end);
    
    for rc2 = rc+1:length(ROIs)
        %%        
        puff2 = ROIs(rc2).puffind;
        signal2 = ROIs(rc2).avgsignal(opt.startPA:end);
        
        if (false)
         plot(tlist,allsignal(rc,:),'b',tlist(puff),signal(puff),'bo')
         hold all
         plot(tlist,allsignal(rc2,:),'r',tlist(puff2),signal2(puff2),'ro')
         hold off
        end
%         
        alldist = pdist2(puff',puff2');
        signaldist(rc,rc2) = min(alldist(:));
        signaldist(rc2,rc) = signaldist(rc,rc2);
    end
end

signaldist(1,1) = 0;

%% get connected clusters of time-similar regions

allpuffmask = zeros(size(ROIs(1).mask));
for rc = 1:length(ROIs)
    if (~isempty(ROIs(rc).puffind))
        allpuffmask = allpuffmask | ROIs(rc).mask;        
    end    
end


% connectivity graph for regions within a cutoff timeseparation between
% closest peak
closereg = signaldist<opt.timesepcut;

congraph = graph(closereg);
regbins = conncomp(congraph);
nbins = max(regbins);

clear peakregind peakregmask peakregcent
ct = 0;
for bc =1:nbins
    reg = find(regbins==bc);
    
    regmask = zeros(size(ROIs(1).mask));
    for rc = reg
        regmask = regmask | ROIs(rc).mask;
    end
    
    %imshowpair(CL.ERimg,regmask)
    
    % break up into spatial connected components
    spcomp = bwconncomp(regmask>0);
    
    % mask for the connected region corresponding to these time peaks
    for sc = 1:spcomp.NumObjects;
        ct =ct+1;
        ids = spcomp.PixelIdxList{sc};
        compmask = zeros(size(allpuffmask));
        compmask(ids) = 1;
        
        peakregmask(:,:,ct) = compmask;
        % index among compared signal bins
        peakregind(ct) = bc;
        
        % define centroid
        tmp = regionprops(peakregmask(:,:,ct),'Centroid');
        cent = tmp.Centroid;
        peakregcent(ct,:) = cent;
    end
end
    

