function [peakregcent,peakregmask] = spatialPeakReg(peakregmask,imgs,peakoptions)
% cluster together ROIs with similar peak positions in time
% that are also connected in space
% peakoptions = options for peak finding

% returns new center points for the peak regions
% corresponding to the highest temporal peak in that region

opt = struct();
opt.dodisplay = 0;
opt.puffsize = 10;
opt.Imax = 0.0025;
opt.ppmass = 95;
opt.roiwidth=20;

if (exist('options','var'))
    opt = copyStruct(options,opt,1);
end

dodisplay=opt.dodisplay;
makerois = false;

nreg = size(peakregmask,3);

% get traces for each region
imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
maskshape = reshape(peakregmask,size(imgs,1)*size(imgs,2),nreg);
regionTraces = (maskshape'*imgshape);

%[keepreg, newROIs] = filterForPeaks(newROIs,options);

% structures representing the rois for the regions (mask and signal only)
clear regROIs
for rc = 1:nreg
    regROIs(rc) = struct('avgsignal',regionTraces(rc,:), 'mask',peakregmask(:,:,rc));
end
[keepreg, regROIs] = filterForPeaks(regROIs,peakoptions);

regROIs = regROIs(keepreg);
peakregmask = peakregmask(:,:,keepreg);


nreg = length(regROIs);
%%
if (dodisplay)
    cmat = jet(nreg)       
    for rc = 1:nreg
        puffind = regROIs(rc).puffind;
        plot(regROIs(rc).avgsignal,'Color',cmat(rc,:))        
        hold all
        plot(regROIs(rc).basesignal,'Color',cmat(rc,:))
        plot(puffind,regROIs(rc).avgsignal(puffind),'o','Color',cmat(rc,:))
    end
    hold off
end
%%

for rc = 1:nreg
    %% specific region mask
    mask = regROIs(rc).mask;
    % dilate max
    se = strel('disk',round(opt.puffsize/2));
    dilmask = imdilate(mask,se);

    % get time range for largest peak (gaussian filter)
    [h,maxind] = max(regROIs(rc).puffheights);
    puffind = regROIs(rc).puffind(maxind);

    %% time smooth images around peak

    % make a gaussian window filter
    w = ceil(regROIs(rc).puffwidths(maxind)/2);% half-width window
    gw = gausswin(w*2+1,w);
    % get weighted time average 
    inds = puffind(1)-w:puffind(1)+w;
    permimg= permute(imgs(:,:,inds),[3,1,2]);
    imgtimefilt = squeeze(pagemtimes(gw',permimg));

    %% spatially smooth image with averaging filter
    Kaverage = filter2(fspecial('average',opt.puffsize),imgtimefilt)/255;
    imgsmoothed = Kaverage.*dilmask;   

    % parameters for feature finding
    Imax= opt.Imax; % percentile cutoff for peaks
    massmin = prctile(imgsmoothed(:),opt.ppmass)*(2*opt.puffsize)^2;

    [goodpts,tmp,pts] = detectParticles(imgsmoothed,opt.puffsize,'Imax',Imax,'massmin',massmin,'massmax',inf,...
      'maxecc',1);
    % keep only those points within the original mask
    maskvals = interp2(mask,goodpts(:,1),goodpts(:,2));
    goodpts = goodpts(abs(maskvals-1)<eps,:);
    goodpts = sortrows(goodpts,3,'descend');
    peakcent = goodpts(1,1:2); % save the eak with the biggest mass


    if (dodisplay)
        imshowpair(imgtimefilt,imgsmoothed)
        hold all
        plot(goodpts(:,1),goodpts(:,2),'y.','MarkerSize',10)
        hold off       
    end

    peakregcent(rc,:) = peakcent;

    if (makerois)
        %% make an ROI around the peak
        rect = [peakcent-[opt.roiwidth/2,opt.roiwidth/2] opt.roiwidth opt.roiwidth];
        if (dodisplay)
            imshow(imgtimefilt,[])
            hold all
            plot(peakcent(1),peakcent(2),'m.','MarkerSize',10)
            peakroi = drawrectangle('Position',rect,'Color','g')
            peakroimask = createMask(peakroi);
            hold off
        end

        %% plot signal for new ROI
        if (dodisplay)

            % reshape each image into a single vector
            imgshape = reshape(imgs,size(imgs,1)*size(imgs,2),size(imgs,3));
            maskshape = reshape(peakroimask,size(imgs,1)*size(imgs,2),size(peakroimask,3));

            regionTracePeak = (maskshape'*imgshape);
            tvals = 1:length(regionTracePeak);

            plot(tvals,regionTracePeak)
            hold all
            plot(tvals(puffind(1)),regionTracePeak(puffind(1)),'ro')
            hold off
        end
    end
end