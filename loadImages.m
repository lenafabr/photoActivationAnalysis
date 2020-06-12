function [imgs,startPA,rawimgs] = loadImages(dirname,prefile,PAfile,options)
% load in image series from tif stacks
% --------
% Input:
% ---------
% dirname = directory name (applied as prefix to both file names)
% prefile has filename for images before photoactivation 
% PAfile has filename for images once photoactivation starts
% can pass '' for either prefile or PAfile, to only load one stack
% -----------
% options = structure with extra options (ignore if just using defaults)
% options.denoise = filter to use for denoising (default = 'median'; pass
% '' if you want no denoising)
% options.medpx = pixels to use for median filter (default=3)
% -------
% Output:
% -------
% imgs = 3d matrix, imgs(:,:,fc) has the image for the fc-th frame
% images are stacked with the ones from prefile first, then the ones from
% PAfile
% startPA = index of first image from PAfile
% rawimgs = all the images without denoising

% ----- set default options ----
opt = struct();

opt.denoise = 'median'; % denoising filter
opt.medpx = 3; % pixels to use in median filter

if (exist('options','var'))
    % copy over passed options
    opt = copyStruct(options,opt);
end
%% load in images from prefile
if (isempty(prefile))
    startPA = 1;
    preimgs = [];
else
    filename = [dirname prefile];
    if (~exist(filename,'file'))
        error(sprintf('Pre file not found: %s', filename))
    end
    info = imfinfo(filename); % general info about images
    nimg = length(info);
    presize = [info(1).Height,info(1).Width];
    
    preimgs = zeros(presize(1),presize(2),nimg);
    for fc = 1:nimg
        % load one frame at a time
        im = imread(filename,fc);
        % rescale image values to be btwn 0 and 1
        im = im2double(im);
        preimgs(:,:,fc) = im;
    end
                    
    % index where photoactivated images start
    startPA = nimg + 1;
end 

if (isempty(PAfile))
    PAimgs = [];
else
    filename = [dirname PAfile];
    if (~exist(filename,'file'))
        error(sprintf('Photoactivated file not found: %s', filename))
    end
    info = imfinfo(filename); % general info about images
    nimg = length(info);
    PAsize = [info(1).Height,info(1).Width];
    
    PAimgs = zeros(PAsize(1),PAsize(2),nimg);
    for fc = 1:nimg
        % load one frame at a time
        im = imread(filename,fc);
        % rescale image values to be btwn 0 and 1
        im = im2double(im);
        PAimgs(:,:,fc) = im;
    end
end
    
% concatenate all images
rawimgs = cat(3,preimgs,PAimgs);

%% denoise if desired
imgs = rawimgs;

nimg = size(imgs,3);
for fc = 1:nimg
    % denoise with a median filter over 3 pixels (default)
    if (strcmp(opt.denoise,'median'))
        imgs(:,:,fc) = medfilt2(rawimgs(:,:,fc),[opt.medpx opt.medpx]);
    elseif (~isempty(opt.denoise))
        error(sprintf('Currently the only denoising filters set up are median or nothing. You passed opt.denoise = %s', opt.denoise))
    end
end

end