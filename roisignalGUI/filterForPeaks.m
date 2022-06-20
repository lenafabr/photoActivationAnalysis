function [keepreg,ROIs] = filterForPeaks(ROIs,options)
%% find ROIs that have peaks
% define peaks: peak in detrended signal *and* sufficient negative curvature in smoothed signal
% takes in list of cell ROI structures
% returns: 
% keepreg = list of regions to keep because they have puffs
% ROIs = roi structures updated to contain info on puffs

% default options
opt = struct();


% polynomial fit order for detrending
opt.detrendorder = 10;

% svg order and window size for smoothing signal to get curvature
opt.order = 3;
opt.framelen = 15;
% starting frame for activated signal
opt.startPA = 1;

% scaling factor for universal threshold
% for defining peak in raw data
opt.cutscalepeak = 1;
% for defining sufficient negative curvature in smoothed data
opt.cutscalecurvature = 1;

if (exist('options','var'))
    % copy over passed options
    opt = copyStruct(options,opt);
end

% ---------------
[b,g] = sgolay(opt.order,opt.framelen);

dodisplay = 1;
prevmax = 0;
keepreg =[];
cmap = jet(length(ROIs));
clear dervsig allsignal maxpuffheight maxpuffwidth maxpuffderv

warning('off','signal:findpeaks:largeMinPeakHeight')

xvals0= linspace(-5,5,100);
gauss0 = normpdf(xvals0)*sqrt(2*pi);
for rc = 1:length(ROIs)
    signal = ROIs(rc).avgsignal(opt.startPA:end);
    nt = length(signal);
    tlist = 1:length(signal);
    
    % smooth signal, for detrending
    [p,s,mu] = polyfit((1:length(signal)),signal,opt.detrendorder);
    f_y = polyval(p,(1:numel(signal))',[],mu);
    %f_y = smoothdata(signal,'sgolay',length(signal)/10)';
    
    % find peaks
    detrendsignal = signal-f_y';    

    cutoff = opt.cutscalepeak*mad(detrendsignal,1)/0.6745*sqrt(2*log(length(detrendsignal)));    
    
    [~,puffind,puffw] = findpeaks(detrendsignal,'MinPeakProminence',cutoff,'MinPeakHeight',cutoff);    
    
    % estimate smoothed signal and first 2 derivs
    dx = zeros(length(detrendsignal),3);
    dt = 1; % time step
    for p = 0:2
        dervsig(:,p+1) = conv(detrendsignal, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end

    smcutoff = opt.cutscalecurvature*mad(dervsig(:,3),1)/0.6745*sqrt(2*log(size(dervsig,1) - opt.framelen));  
    goodind = find(dervsig(puffind,3)<-smcutoff);
    goodpuff = puffind(dervsig(puffind,3)<-smcutoff);

    % check for low-signal region
    %if (islowsig(rc)); goodpuff = []; end
    
    ROIs(rc).puffind = goodpuff;
    %allpuffreg{rc} = rc*ones(size(goodpuff)); % which region these puffs are in
    %allsignal(rc,:) = signal;        

    if (~isempty(goodpuff))                
        %shift = (prevmax-prctile(signal,10));
        %plot(tlist,signal+shift,tlist,f_y+shift,tlist(goodpuff),signal(goodpuff)+shift,'ro','Color',cmap(rc,:))        
        %title(sprintf('ROI %d',rc))
        %prevmax = max(signal+shift);
        %hold all
        %text(tlist(end),signal(end)+shift,sprintf('%d',rc),'Color',cmap(rc,:))
        
        ROIs(rc).maxpuffheight = max(signal(goodpuff)/cutoff);
        ROIs(rc).maxpuffderv = max(-dervsig(goodpuff,3)/smcutoff);
        ROIs(rc).maxpuffwidth = max(puffw(goodind));
        
        ROIs(rc).puffwidths = puffw(goodind);
        
        % get puffintegrals
        peakinteg = zeros(1,length(goodpuff));
        for cc = 1:length(goodpuff)            
            w = puffw(goodind(cc));
            base = f_y(goodpuff(cc));
            %gauss = gauss0*(signal(goodpuff(cc))-base) + base;
            %xvals = xvals0*w/2 + tlist(goodpuff(cc));        
            %plot(xvals,gauss+shift,'--','Color',cmap(rc,:))
        
            % integral of peak            
            peakinteg(cc) = (signal(goodpuff(cc))-base)*sqrt(2*pi*w/2);
        end
        ROIs(rc).puffinteg = peakinteg; 
        ROIs(rc).basesignal = f_y;
        
        keepreg(end+1) = rc;
    end
end
hold off
%title(opt.Name)