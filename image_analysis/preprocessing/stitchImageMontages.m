function upperleft = stitchImageMontages(dataDir, opts)
%function for stitching image grids from the Nikon or Dragonfly microscope
%stitching can be done while image acquisition is taking place, or all at
%once afterwards

%additional functionality to add:
%   -option to stitch only a subset of the channels?
%   -optionally base the names of stitched files on those of the raw data;
%    for now parse the filenames and save the format in a text file for
%    reference***
%   -save upperleft to reconstruct stitching without needing to recompute
%    overlap later?
%   -add dicroic correction -> this was in the code in stitchedPreviews
%   -determine variable class from loaded images -> unnecessary for now
%   because our images are always acquired as 16 bit integers (uint16)


%% parse options

bare = '_p%.4d_w%.4d_t%.4d.tif';
mipdir = fullfile(dataDir, 'MIP'); %folder in which to save MIPs
if ~exist(mipdir,'dir'), mkdir(mipdir); end

load(fullfile(dataDir,'meta.mat'),'meta')

previewdir = fullfile(dataDir,'previews'); %folder in which to save previews
if ~exist(previewdir,'dir'), mkdir(previewdir); end

if isfield(opts,'ext') && ~isempty(opts.ext)
    %manually specify the image file extension
    ext = opts.ext;
else
    %look for files from a list of microscope image file formats
    exts = {'.ims','.nd2','.lif'};
    idx = 1;
    criterion = true;
    while criterion && idx <= length(exts)
        ext = exts{idx};
        listing = dir(fullfile(dataDir,['*',ext]));
        if ~isempty(listing)
            criterion = false;
        end
        idx = idx + 1;
    end
end

if isfield(opts,'nseconds')
    %amount of time in seconds to wait for updates if stitching is being 
    %done while images are taken and transferred
    nseconds = opts.nseconds;
else
    nseconds = 60;
end

if isfield(opts,'nucChannel')
    nucChannel = opts.nucChannel;
else
    warning('no nuclear channel specified, assuming first channel')
    nucChannel = 0;
end

if isfield(opts,'ngrids')
    ngrids = opts.ngrids;
else
    ngrids = [];
end

if ~isfield(opts,'tmax') || isempty(opts.tmax)
    opts.tmax = Inf;
end

if isfield(opts,'montageFusionGrid')
    %boolean for if the matlab script multiMontageFusion was used to 
    %generate positions for the grids
    montageFusionGrid = opts.montageFusionGrid;
else
    montageFusionGrid = false;
end

if isfield(opts,'mipext')
    mipext = opts.mipext;
else
    %should this be jpg by default, or tiff?
    mipext = '.jpg';
end

if isfield(opts,'stitchmode')
    stitchmode = opts.stitchmode;
else
    stitchmode = 'weightedsample';
end

gridsize = opts.gridsize;
ppg = gridsize(1)*gridsize(2); %positions per grid
montageOverlap = opts.montageOverlap; %percent overlap between adjacent images

%flat field correction
if isfield(opts,'doFlatFieldCorrection')
   doFlatFieldCorrection = opts.doFlatFieldCorrection;
else
    doFlatFieldCorrection = false;
end
if doFlatFieldCorrection
    if ~strcmp(ext, '.nd2')
        doFlatFieldCorrection = false;
        warning('flat field correction can only be used with Nikon files, option set to false')
    else
        %need to add these matrices to the lab repo
        load('correctionMats.mat','Gps','D')
        %should make ffchannels part of the data to be loaded here in case
        %we add more channels (i.e. yellow, cyan) later
        ffchannels = {'405-DGRI','488-DGRI','561-DGRI','640-DGRI'};
    end
end

%%% preview video options %%%
if ~isfield(opts,'makePreviewVideos')
    makePreviewVideos = true;
else
    makePreviewVideos = opts.makePreviewVideos;
end

%option to crop the preview video by a given number of pixels?
if ~isfield(opts, 'margin')
    margin = 0;
else
    margin = opts.margin;
end

% subsampling factor -> image is reduced in size
if ~isfield(opts, 'ss')
    if meta.xSize <= 1024
        ss = 2;
    else
        ss = 4;
    end
else
    ss = opts.ss;
end
margin = round(margin/ss);

if ~isfield(opts,'previewFormat')
    opts.previewFormat = 'mp4';
end

if ~isfield(opts,'FrameRate')
    opts.FrameRate = 5;
end

%% automatic settings
%if the number of grids/colonies is empty, determine it automatically
if isempty(ngrids)
    if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
        listing = dir(fullfile(dataDir,['*',ext]));
        r = bfGetReader(fullfile(dataDir,listing(1).name));
        %do we also need to account for multiple files with multiple series?
        npos = r.getSeriesCount;
        seriesflag = true;
        
        if npos == 1
            npos = length(listing);
            seriesflag = false;
        end
    else
        npos = length(listing);
        seriesflag = false;
    end
    ngrids = npos/ppg;
else
    seriesflag = false;
end

%determine order of images in a stitched grid
order = NaN(ppg,2);
idx = 1;
if ~montageFusionGrid
    for ii = 1:gridsize(1)
        if mod(ii,2) == 0 %alternate direction for odd and even rows
            roworder = gridsize(2):-1:1;
        else
            roworder = 1:gridsize(2);
        end
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
else
    for ii = 1:gridsize(1)
        roworder = 1:gridsize(2);
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
end

grid = NaN(gridsize);
for idx = 1:ppg
    ii = order(idx,1); jj = order(idx,2);
    grid(ii,jj) = idx;
end
disp('grid layout')
disp(grid)

fprintf("number of stitched grids = %d\n",ngrids)


%% do stitching
%iterate over grids; structured as a while loop instead of a for loop
%iterating over grids for the case in which we are doing stitching as the
%images are acquired and transferred
pidx = 1;
while pidx <= ngrids %TODO: add the option to only do stitching for specified colonies/grids
    %find all the image (nd2) files in the directory
    listing = dir(fullfile(dataDir,['*',ext]));
    nfiles = numel(listing);
    if ~seriesflag
        fprintf('%d of %d files found\n',nfiles,ngrids*ppg)
    end
    %if there are at least enough images to stitch the pidxth grid, stitch
    %that one and increment pidx
    if nfiles >= pidx*ppg || seriesflag
        if pidx == 1
            %assume the number of time points is the same across all files
            %(for this type of data in general it should be = 1)
            r = bfGetReader(fullfile(dataDir,listing(1).name));
            nt = r.getSizeT;
            nt = min(nt,opts.tmax);
        end
        fprintf('Stitching grid #%d of %d\n',pidx,ngrids)
        
        for ti = 1:nt
            fprintf('t = %d of %d\n',ti,nt)
            imgs = cell(gridsize);
            nucmips = cell(gridsize);
            for idx = 1:ppg
                %overall index of the (idx)th image in the (pidx)th grid
                fidx = (pidx - 1)*ppg + idx;
                fprintf('.')
                
                if seriesflag
                    %if all positions are saved as series in a single file
                    r.setSeries(fidx-1);
                else
                    fname = fullfile(dataDir,listing(fidx).name);
                    r = bfGetReader(fname);
                end
                
                if idx == 1
                    m = r.getSizeY; n = r.getSizeX; nz = r.getSizeZ; nc = r.getSizeC;
%                     nt = r.getSizeT; %set nt above to allow us to set tmax more easily
                    pixelOverlapY = round(m*montageOverlap/100); pixelOverlapX = round(n*montageOverlap/100);
                    stitchedSize = round([gridsize(1) - montageOverlap*(gridsize(1)-1)/100,...
                        gridsize(2) - montageOverlap*(gridsize(2)-1)/100].*[m,n]);
                    omeMeta = r.getMetadataStore();
                    channelNames = cell(1,nc);
                    for ci = 1:nc 
                        channelNames{ci} = char(omeMeta.getChannelName(0,ci-1));
                    end
                end
                
                
                ii = order(idx,1); jj = order(idx,2);
                img = zeros(m,n,nz,nc,'uint16');
                for ci = 1:nc
                    for zi = 1:nz
                        im = bfGetPlane(r, r.getIndex(zi-1,ci-1,ti-1)+1);
                        img(:,:,zi,ci) = im;
                    end
                end
                imgs{ii,jj} = img;
                nucmips{ii,jj} = max(squeeze(img(:,:,:,nucChannel+1)),[],3);
            end
            fprintf('\n')

            if ti == 1 %double check the overlap values in pixels on the first time point
                %check pixelOverlapY
                [~,I] = max((cellfun(@(x) mean(x(end-pixelOverlapY+1:end,:),'all'),nucmips(1:end-1,:))),[],'all','linear');
                [ii,jj] = ind2sub(gridsize - [1 0],I);
                im1 = nucmips{ii,jj}(end-pixelOverlapY+1:end,:);
                im2 = nucmips{ii+1,jj}(1:pixelOverlapY,:);
                [offset, ~, ~, ~] = xcorr2fft(im1, im2);

                if abs(offset/pixelOverlapY) > 0.25
                    warning('detected pixel overlap very different from nominal value, sticking with nominal value');
                else
                    pixelOverlapY = pixelOverlapY + offset;
                end

                %check pixelOverlapX
                [~,I] = max((cellfun(@(x) mean(x(:,end-pixelOverlapX+1:end),'all'),nucmips(:,1:end-1))),[],'all','linear');
                [ii,jj] = ind2sub(gridsize - [0 1],I);
                im1 = nucmips{ii,jj}(:,end-pixelOverlapX+1:end);
                im2 = nucmips{ii,jj+1}(:,1:pixelOverlapX);
                [~, offset, ~, ~] = xcorr2fft(im1, im2);
                
                if abs(offset/pixelOverlapX) > 0.25
                    warning('detected pixel overlap very different from nominal value, sticking with nominal value');
                else
                    pixelOverlapX = pixelOverlapX + offset;
                end
                
                %initialize a stack of subsampled images for making preview videos
                if nt > 1 && makePreviewVideos
                    smallsize = ceil(stitchedSize/ss);
                    preview = zeros([smallsize,nc,nt],'uint16');
                end
            end

            upperleft = registerImageGrid_v3(nucmips, [pixelOverlapY pixelOverlapX]);
            if strcmp(stitchmode,'weightedsample')
                zidxs = randomIndicesDistWeighted(upperleft, nucmips);
            end
            
            for ci = 1:nc
                suffix = sprintf(bare,pidx-1,ci-1,ti-1);
                fname = fullfile(dataDir,['stitched',suffix]);
                mipname = fullfile(mipdir,['stitched_MIP',suffix]);
                mipname = strrep(mipname,'.tif',mipext); %save the MIP as a png
                previewname = fullfile(previewdir,['stitched_preview',suffix(1:end-4),'.jpg']);

                stitched = zeros(stitchedSize(1),stitchedSize(2),nz,'uint16');
                for zi = 1:nz
                    if zi == 1
                        mode = 'overwrite';
                    else
                        mode = 'append';
                    end
                    ims = cellfun(@(x) x(:,:,zi,ci),imgs,'UniformOutput',false);
                    
                    if doFlatFieldCorrection
                        if ci == 1 && zi == 1 && ti == 1
                            disp('using flat field correction')
                        end
                        I = find(strcmp(ffchannels,channelNames{ci}));
                        if ~isempty(I)
                            G = Gps{I};
                            ims = flatFieldCorrection(ims,G,D);
                        else
                            warning('cannot do flat field correction for the current channel')
                            disp(strcat("channel name = ",channelNames{ci}))
                        end
                    end
                    
                    if strcmp(stitchmode,'distweight')
                        [stitch, ~] = stitchImageGridDistWeighted(upperleft, ims);
                    elseif strcmp(stitchmode,'intensityweight')
                        [stitch, ~] = stitchImageGridIntensityWeighted(upperleft, ims);
                    elseif strcmp(stitchmode,'randomsample')
                        [stitch, ~] = stitchImageGridRandSample(upperleft, ims);
                    elseif strcmp(stitchmode,'weightedsample')
                        stitch = applyRandomIndices(upperleft, ims, zidxs);
                    end
                    
                    %trim to projected size
                    tmp = zeros(stitchedSize,'uint16');
                    yrange = 1:min(stitchedSize(1),size(stitch,1));
                    xrange = 1:min(stitchedSize(2),size(stitch,2));
                    tmp(yrange, xrange) = stitch(yrange, xrange);
                    stitched(:,:,zi) = tmp;

                    imwrite(tmp,fname,'WriteMode',mode)
                end
                %write mip
                [MIP,~] = max(stitched,[],3);
                if strcmp(mipext,'.jpg')
                    imwrite(im2double(MIP),mipname,'Quality',99)
                else
                    imwrite(MIP,mipname)
                end
                
                %write preview image
                if ti == 1 || ti == nt
                    imwrite(im2double(imadjust(MIP,stitchedlim(MIP))),previewname)
                end
                %subsample for preview video
                if nt > 1 && makePreviewVideos
                    small = imfilter(MIP,ones(ss)/ss^2); %average pixel values in an ss x ss box
                    small = small(1:ss:end,1:ss:end); %subsample averaged values
                    preview(:, :, ci, ti) = small;
                end
            end
        end
        
        %%% write preview video in each channel %%%
        if nt > 1 && makePreviewVideos
            fprintf('writing preview videos for grid %d\n',pidx)
            for ci = 1:nc
                % set lookup table -> find the average intensity of the
                % image in this channel at each time, then use the image in
                % the top 98th percentile to adjust the contrast
                timepointavgs = squeeze(mean(preview(:,:,ci,:),[1,2]));
                [~,I] = sort(timepointavgs,'descend');
                Ilim = stitchedlim(preview(:,:,ci,I(ceil(nt/50))));

                % crop
                preview = preview(1+margin:end-margin, 1+margin:end-margin,:,:);
                
                %filename of the saved video
                videoname = fullfile(previewdir,['stitchedPreview',sprintf('_p%.4d_w%.4d',pidx-1,ci-1),'.',opts.previewFormat]);
                
                % save as compressed video
                if strcmp(opts.previewFormat,'avi')
                    v = VideoWriter(videoname,'Uncompressed AVI'); %#ok<TNMLP>
                elseif strcmp(opts.previewFormat,'mp4')
                    v = VideoWriter(videoname,'MPEG-4'); %#ok<TNMLP>
                else
                    error('unknown video format');
                end
                
                v.FrameRate = opts.FrameRate;
                v.Quality = 100;
                open(v)
                for ti = 1:nt
                    writeVideo(v,im2double(imadjust(preview(:,:,ci,ti),Ilim)))
                end
                close(v);
            end
        end
        
        %next grid
        pidx = pidx + 1;
    end
    %pause if there are not already enough files to stitch the next grid
    if numel(listing) < ngrids*ppg && numel(listing) < pidx*ppg && ~seriesflag
        disp('paused')
        timestamp=datetime('now');
        disp(timestamp)
        pause(nseconds)
    end
end

%update and save metadata
ncond = meta.nWells;
meta.posPerCondition = ngrids/ncond;
meta.nPositions = ngrids;
meta.montageGridSize = gridsize;
meta.montageOverlap = montageOverlap;
meta.xSize = stitchedSize(2);
meta.ySize = stitchedSize(1);

save(fullfile(dataDir,'meta.mat'),'meta');





end


















