classdef Position < handle
    % Data class to store cell data in a field of view

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        ncells              % number of cells

        cellData            % structure array indexed by time:
                            % - XY
                            % - area
                            % - nucLevel:       nCell x nChannels array
                            % - cytLevel
                            % - background:     nChannels long vector

        timeTraces          % structure: reorganizes cellData
                            % - nucLevelAvg:    vector indexed by time
                            % - cytLevelAvg
                            % - background
                            % NOT IMPLEMENTED YET: traces of tracked cells

        dataChannels        % we may not want to quantify all channels

        nChannels           % number of colors imaged
        nucChannel          % nuclear channel
        
        nTime               % number of time points
        filenameFormat
        tree                %stores a tree object with tracks (deprecated)
        G                   %digraph with links from tracking
    end

    properties (SetAccess = protected)
        ID                  % identifier of colony
    end

    properties (Dependent)

        filename
        barefilename;       % filename without extension
        %data                % cell by cell data for colony
                            % data cols: x, y, area, -1,
                            % (nuclear, cytoplasmic) mean for each channel
    end

    methods

        % constructor
        function this = Position(varargin)
            % Position(nChannels, filenameFormat, nTime, nucChannel)
            % Position(meta)
            % Position(meta, id)

            if nargin == 4
                nChannels = varargin{1};
                filenameFormat = varargin{2};
                nTime = varargin{3};
                nucChannel = varargin{4};

            elseif nargin == 1 || nargin == 2
                
                meta = varargin{1};
                nChannels = meta.nChannels;
                filenameFormat = meta.filenameFormat;
                nTime = meta.nTime;
                nucChannel = meta.nucChannel;
                
                if nargin == 2
                    this.ID = varargin{2};
                    if numel(meta.posPerCondition) > 1
                        conditionIdx = find(this.ID >= meta.conditionStartPos,1,'last');
                        conditionPos = this.ID - meta.conditionStartPos(conditionIdx) + 1;
                    else % backwards compatibility
                        conditionIdx = ceil(this.ID/meta.posPerCondition);
                        conditionPos = rem(this.ID-1, meta.posPerCondition)+1;
                    end
                    if isempty(filenameFormat) && ~isempty(meta.fileNames{1})
                        filenameFormat = meta.fileNames{conditionIdx}{conditionPos};
                    end
                end

            elseif nargin == 0 % matlab sucks, always call superclass constructor with no args
                return
            else
                error('wrong number of arguments');
            end

            if ~exist('nTime','var')
                nTime = 1;
                this.timeTraces = [];
            else
                this.timeTraces = struct();
            end

            if ~isnumeric(nChannels)
                error('first argument is nChannels, which is a number');
            end
            if ~isnumeric(nTime)
                error('third argument is nTime, which is a number');
            end

            this.nucChannel = nucChannel;
            this.nChannels = nChannels;
            this.nTime = nTime;

            % strip off path in case it was provided
            if ~isempty(filenameFormat)
                [~, name, ext] = fileparts(filenameFormat);
                this.filenameFormat = char(strcat(name,ext));
            end

            this.cellData = struct();
        end

        % saving and loading
        %---------------------------------

        function img = loadImage(this, dataDir, channels, time)
            % load position image
            %
            % img = loadImage(dataDir, channels)
            % img = loadImage(dataDir, channels, time)
            %
            % dataDir:  main data directory
            % channels: desired channels to be loaded, leave empty for all
            % time:     time to load
            %
            % img:      loaded image

            if ~exist('channels','var')
                channels = 1:this.nChannels;
            end
            if ~exist('time','var')
                time = 1;
            end
            if ~exist('useMIP','var')
                useMIP = false;
            end

            [~,~,ext] = fileparts(this.filenameFormat);

            if strcmp(ext,'.tif') || strcmp(ext,'.btf')

                s = strsplit(this.filenameFormat,'_F[0-9]','DelimiterType','RegularExpression');
                if numel(s) > 1
                    fnameFormat = [s{1} '_MIP_p%.4d_w%.4d.tif'];
                else
                    fnameFormat = this.filenameFormat;
                end

                for cii = 1:numel(channels)

                    if contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
                        fname = sprintf(fnameFormat,this.ID - 1, channels(cii), time-1);
                        iminf = imfinfo(fullfile(dataDir,fname));
                        im = zeros([iminf(1).Height iminf(1).Width numel(iminf)],'uint16');
                        for zii = 1:numel(iminf)
                            im(:,:,zii) = imread(fullfile(dataDir, fname), zii);
                        end
                        img{cii} = im;
                        
                    elseif ~contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
                        % this doesn't handle z and t in the same tiff
                        fname = sprintf(fnameFormat,this.ID - 1, channels(cii));
                        img{cii} = imread(fullfile(dataDir, fname), time);
                        
                    else
                        % this doesn't handle z and t in the same tiff
                        fname = sprintf(fnameFormat,this.ID - 1);
                        img{cii} = imread(fullfile(dataDir, fname), time);
                    end
                end
                img = cat(3,img{:});

            % bioformats
            elseif strcmp(ext,'.vsi') || strcmp(ext, '.oif') ...
                    || strcmp(ext, '.oib') || strcmp(ext, '.ims') || strcmp(ext, '.nd2')

                fname = fullfile(dataDir, this.filename);

                r = bfGetReader(fname);
                img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ()], 'uint16');
                for cii = 1:numel(channels)
                    for zi = 1:r.getSizeZ()
                        img(:,:,cii,zi) = bfGetPlane(r, r.getIndex(zi-1,channels(cii),time-1)+1);
                    end
                end
                r.close();

                img = squeeze(img);
            else
                error('unknown format');
            end
            if time == 1
                disp(['loaded image ' fname]);
            end
        end

        function seg = loadSegmentation(this, dataDir, channel, segfileformat)
            % load segmentation
            %
            % seg = loadSegmentation(dataDir, channel)
            %
            % dataDir:  main data directory
            %
            % seg:      binary image of nuclei
            %
            % segmentation is expected to have same name as raw data
            % but with .tif -> '_Simple Segmentation.h5'
            %
            % by convention use green for foreground in Ilastik

            if ~exist('channel','var')
                error('please provide channel');
            end

            if exist('segfileformat','var')

                fname = sprintf(segfileformat, channel);
            else

                % guessing filename based on conventions
                %------------------------------------------

                % possible cases for filename vs segmentation:
                % 1) 1.a) Andor multiposition
                %    blabla_F00.ims -> blabla_MIP_p0000_w0000.tif 
                %                   -> blabla_MIP_p0000_w0000_Simple Segmentation.h5
                %    1.b) can use MIP instead of raw data for quantification
                %
                % 2) Snaps with well names for condition
                %    A1_blabla.vsi  -> A1_blabla_MIP_w0000.tif
                %                   -> A1_blabla_MIP_w0000_Simple Segmentation.h5
                
                barefname = {};

                % case 1:
                
                % split by _ <letter> <number>
                s = strsplit(this.filename,'_[Ffwptm][0-9]+','DelimiterType','RegularExpression');
                
                % get the position index from the raw data
                [si,ei] = regexp(this.filename,'_[Ffpm][0-9]+');
                pi = str2num(this.filename(si+2:ei));

                % add MIP if the 'raw' data filename is not a MIP (usual
                % case) because we always segment MIPS
                % 'raw' data is MIP in the case for the stitched grids
                if contains(s{1},'_MIP') % case 1.b)
                    barefname{1} = sprintf([s{1}  '_p%.4d'], pi);
                else  % case 1.a)
                    barefname{1} = sprintf([s{1}  '_MIP_p%.4d'], pi);
                end
                
                % case 2:
                barefname{2} = this.filename(1:end-4);
                
                % find the first barefname for which there is an h5 file
                i = 0;
                listing = [];
                while(isempty(listing))
                    i = i+1;
                    if i > numel(barefname)
                        error(['segmentation for ' this.filename ' for channel ' num2str(channel) ' not found in ' dataDir]);
                    else
                        listing = dir(fullfile(dataDir,[barefname{i} '*.h5']));
                    end
                end

                fname = [];
                for i = 1:numel(listing)
                    % DON'T COMMENT THE IF STATEMENT BELOW, IT WILL LOAD THE WRONG
                    % CHANNEL
                    if ~isempty(strfind(listing(i).name, sprintf('_w%.4d',channel)))...
                            || ~isempty(strfind(listing(i).name,sprintf('_c%d',channel)))

                        fname = fullfile(dataDir,listing(i).name);
                    end
                end
            end

            % actually reading the file
            %----------------------------

            if ~exist(fname,'file')

                warning(['segmentation for channel ' num2str(channel) ' not found in ' dataDir, ', may be naming convention problem']);
                seg = [];

            else
                seg = h5read(fname, '/exported_data');
                
                % Probabilities -> treshold to get binary
                % this is not set up to deal with more than 2 classes yet
                % (fg/bg)
                if size(seg,1) > 1
                    ssegsize = size(seg);
                    ssegsize(1) = 1;
                    simpleseg = false(ssegsize);
                    simpleseg(seg(1,:,:,:) >= 0.5) = true;
                    seg = squeeze(simpleseg);
                    
                % Simple Segmentation (could be more than 2 classes)
                else
                    seg = squeeze(seg);
                end
                % Ilastik output data has xy transposed
                seg = permute(seg, [2 1 3]);
                disp(['loaded segmentation ' fname]);
            end
        end

        function MIPidx = loadMIPidx(this, dataDir, channel, time)
            % load MIP index of position
            %
            % MIPidx = loadMIPidx(this, dataDir, channel)
            %
            % dataDir:  main data directory
            % channels: desired channels to be loaded, leave empty for all
            % time:     time to load, default 1 for static images
            %
            % img:      loaded image

            if ~exist('time','var')
                time = 1;
            end

            fname = this.filename;
            startIndex = regexp(fname,'_[Ffwptm][0-9]+');
            barefname = fname(1:startIndex(1));

            listing = dir(fullfile(dataDir,[barefname '*MIPidx*']));
            if contains(listing(1).name, '_t')
                fname = sprintf([barefname 'MIPidx_p%.4d_w%.4d_t%.4d.tif'], this.ID-1, channel, time-1);
            else
                fname = sprintf([barefname 'MIPidx_p%.4d_w%.4d.tif'], this.ID-1, channel);
            end
            fname = fullfile(dataDir, fname);

            % when I make MIPs I replace m or f by p, I want all positions
            % labeled the same and thats how Andor should have done it
            montageIdx = regexp(fname,'(m|f)[0-9]+','once');
            if ~isempty(montageIdx)
                fname(montageIdx) = 'p';
            end
            
            if ~exist(fname, 'file')

                error(['MIPidx file does not exist: ' fname]);
            else
                if contains(listing(1).name, '_t')
                    MIPidx = imread(fname);
                else
                    MIPidx = imread(fname,time);
                end
            end

            %disp(['loaded MIPidx ' fname]);
        end

        % process data
        %---------------------------------
        
        function addCellLabels(this, labelMatrix)

            for ti = 1:size(labelMatrix,3)
                if this.ncells(ti) > 0
                    L = labelMatrix(:,:,ti);

                    XY = this.cellData(ti).XY;
                    i = round(XY(:,2));
                    j = round(XY(:,1));
                    lini = sub2ind(size(L),i,j);
                    labels = L(lini);
                    this.cellData(ti).labels = labels;
                else
                    this.cellData(ti).labels = [];
                end
            end
        end

        function debugInfo = extractData(this, dataDir, opts)
            % populate the data array
            %
            % extractData(dataDir, nuclearChannel)
            % extractData(dataDir, nuclearChannel, options)
            %
            % dataDir:              data directory
            %
            % options:              struct with fields
            %
            % -dataChannels         channels for which to extract data
            %                       default is all
            % -tMax                 maximal time, default nTime
            % -------------------------------------------------------------
            % -cleanupOptions       options for nuclearCleanup
            % -dirtyOptions         if using dirtyNuclearSegmentation
            % -nucShrinkage         number of pixels to shrink nuclear mask
            %                       by to get more accurate nuclear
            %                       readouts
            % -cytoplasmicLevels    extract approximate cytoplasmic levels
            % -cytoSize             width of cytoplasmic annulus
            % -cytoMargin           margin around nucleus for cyto ring
            % -fgChannel            channel whose complement is used as a
            %                       background mask
            % -bgMargin             erosion size of the inverse fg
            % -imopenBGsub          open size for imopenBGsub (if using it)
            % -segFG                segmentation foreground channel
            % -------------------------------------------------------------
            % -MIPidxDir            directory of MIPidx files if desired
            % -segmentationDir      directory of segmentation, default
            %                       dataDir
            % -------------------------------------------------------------
            % -nuclearSegmentation  binary stack with 3rd dim time for
            %                       providing a segmentation by hand

            if nargin < 3
                opts = struct();
            end
            if ~isfield(opts,'dataChannels')
                opts.dataChannels = 1:this.nChannels;
            end
            if ~isfield(opts, 'tMax') || opts.tMax > this.nTime
                opts.tMax = this.nTime;
            end
            if ~isfield(opts, 'NCRcutoff') || isempty(opts.NCRcutoff)
                opts.NCRcutoff = Inf*ones([1 numel(opts.dataChannels)]);
            end
            % ----------------------------------------
            if ~isfield(opts, 'cleanupOptions')
                opts.cleanupOptions = struct('separateFused', true,...
                    'clearBorder',true, 'minAreaStd', 1, 'minSolidity',0);
            end
            if ~isfield(opts, 'nucShrinkage')
                opts.nucShrinkage = 0;
            end
            if ~isfield(opts,'cytoplasmicLevels')
                opts.cytoplasmicLevels = false;
            end
            if ~isfield(opts, 'cytoSize')
                opts.cytoSize = 10;
            end
            if ~isfield(opts, 'cytoMargin')
                opts.cytoMargin = 0;
            end
            if ~isfield(opts, 'bgMargin')
                opts.bgMargin = 2;
            end
            if ~isfield(opts, 'imopenBGsub')
                opts.imopenBGsub = [];
            end
            % ----------------------------------------
            if ~isfield(opts, 'segmentationDir')
                opts.segmentationDir = dataDir;
            end
            if ~isfield(opts, 'MIPidxDir')
                opts.MIPidxDir = [];
            end
            if contains(this.filenameFormat,'_MIP') && ~isempty(opts.MIPidxDir)
                error('cannot use MIP instead of z-stack data and use MIPidx at the same time, remove MIPidxDir from options');
            end
            if isfield(opts,'keepPixelIdxs')
                keepPixelIdxs = opts.keepPixelIdxs;
            else
                keepPixelIdxs = false;
            end

            this.dataChannels = opts.dataChannels;
            allChannels = unique([this.nucChannel opts.dataChannels]);

            % for the purpose of the current background subtraction
            if isfield(opts, 'fgChannel')
                allChannels = unique([allChannels opts.fgChannel]);
            else
                bgmask = [];
                fgmask = [];
            end

            % load all available segmentations except for nucleus
            seg = {};
            for ci = 1:numel(allChannels)
                if allChannels(ci) ~= this.nucChannel
                    seg{ci} = this.loadSegmentation(opts.segmentationDir, allChannels(ci));
                    % binarize if only two classes
                    if isfield(opts, 'segFG')
                        seg{ci} = seg{ci} == opts.segFG;
                    end
                end
            end
            % pass nuclear segmentation manually or load Ilastik standard
            nucChannelIdx = find(allChannels == this.nucChannel);
            if isfield(opts, 'nuclearSegmentation')
                seg{nucChannelIdx} = opts.nuclearSegmentation;
                junkmask = [];
                if size(seg{nucChannelIdx},3) ~= this.nTime
                    error('slice number of nuclear segmentation doesnt match number of time points');
                end
            else
                seg{nucChannelIdx} = this.loadSegmentation(opts.segmentationDir, this.nucChannel);
                if max(seg{nucChannelIdx}(:,:,end)) == 3
                    junkmask = seg{nucChannelIdx} == 3;
                else
                    junkmask = [];
                end
                seg{nucChannelIdx} = seg{nucChannelIdx} == opts.segFG;
            end
            if isempty(seg{nucChannelIdx})
                error('nuclear segmentation missing');
            end

            % make it more efficient to deal with e.g. oib
            [~,~,ext] = fileparts(this.filename);
            if strcmp(ext,'.oib')
                disp('.oib: loading all data at once for speed');
                imgdata = readStack2(fullfile(dataDir, this.filename), [], opts.tMax);
            end

            for ti = 1:opts.tMax

                % progress indicator
                fprintf('.');
                if mod(ti,60)==0
                    fprintf('\n');
                end

                % get background mask
                %--------------------------

                if isfield(opts, 'fgChannel') && ~isempty(opts.fgChannel) && ~isempty(seg{opts.fgChannel+1})
                    fgidx = opts.fgChannel + 1;
                    fgmask = imclose(seg{fgidx}(:,:,ti),strel('disk',5));
                    bgmask = imerode(~fgmask,strel('disk', opts.bgMargin));
                    bgmask = bgmask > 0;
                    fgmask = imerode(fgmask,strel('disk', opts.bgMargin));
                else
                    if ti == 1
                        warning('not using background mask');
                    end
                    fgmask = [];
                    bgmask = [];
                end

                % make clean nuclear mask
                %--------------------------

                nucmaskraw = seg{nucChannelIdx}(:,:,ti);
                %nucmaskraw(bgmask) = false;
                if ~isfield(opts, 'nuclearSegmentation') && ~isfield(opts, 'dirtyOptions')

                    if ti == 1, disp('using nuclearCleanup'); end
                    nucmask = nuclearCleanup(nucmaskraw, opts.cleanupOptions);

                elseif isfield(opts, 'nuclearSegmentation')

                    nucmask = nucmaskraw;

                elseif isfield(opts, 'dirtyOptions')

                    disp('using dirty nuclear segmentation');
                    opts.dirtyOptions.mask = nucmaskraw;
                    if strcmp(ext,'oib')
                        imc = squeeze(imgdata(:,:,nucChannelIdx,:,ti));
                    else
                        imc = this.loadImage(dataDir, this.nucChannel, ti);
                    end
                    % dirtyNuclearSegmentation expects a 2D image
                    % (for creating watershed seeds)
                    imc = max(imc,[],3);
                    nucmask = dirtyNuclearSegmentation(imc, opts.dirtyOptions);
                end

                % save clean nuclear mask
                if isfield(opts,'saveCleanNucMask') && opts.saveCleanNucMask

%                     segfname = [this.barefilename '_seg.tif'];
%                     segrawfname = [this.barefilename '_segRaw.tif'];
                    segfname = [this.filename, '0000_seg.tif'];
                    segrawfname = [this.filename, '0000_segRaw.tif'];

                    if ti == 1
                        imwrite(nucmask,fullfile(dataDir, segfname),'Compression','none');
                        imwrite(nucmaskraw,fullfile(dataDir, segrawfname),'Compression','none');
                    else
                        imwrite(nucmask,fullfile(dataDir, segfname),'WriteMode','append','Compression','none');
                        imwrite(nucmaskraw,fullfile(dataDir, segrawfname),'WriteMode','append','Compression','none');
                    end
                end

                % make cytoplasmic mask
                %----------------------

                if opts.cytoplasmicLevels

                    if opts.cytoMargin > 0
                        nucmaskmarg = imdilate(nucmask,strel('disk',opts.cytoMargin));
                    else
                        nucmaskmarg = imerode(nucmask,strel('disk',-opts.cytoMargin));
                    end

                    % watershedding inside the dilation
                    dilated = imdilate(nucmask, strel('disk',opts.cytoSize + opts.cytoMargin));
                    basin = bwdist(dilated);
                    basin = imimposemin(basin, nucmask);
                    L = watershed(basin);
                    
                    nucstats = regionprops(L, 'PixelIdxList');
                    cytCC = struct('PixelIdxList', {cat(1,{nucstats.PixelIdxList})});

                    % replacement of exclusions above
                    for cci = 1:numel(cytCC.PixelIdxList)
                        CCPIL = cytCC.PixelIdxList{cci};
                        CCPIL = CCPIL(dilated(CCPIL));    % exclude outside dilated nuclei
                        if opts.cytoMargin ~= 0
                            CCPIL = CCPIL(~nucmaskmarg(CCPIL)); % exclude margin around nuclei
                        else
                            CCPIL = CCPIL(~nucmaskraw(CCPIL));% exclude nuclei before cleanup
                            CCPIL = CCPIL(~nucmask(CCPIL));% exclude nuclei 
                        end
                        if ~isempty(fgmask)
                            CCPIL = CCPIL(fgmask(CCPIL));% exclude background
                        end
                        if ~isempty(junkmask)
                            CCPIL = CCPIL(~junkmask(CCPIL));
                        end
                        cytCC.PixelIdxList{cci} = CCPIL;
                    end

                    % make it a proper CC struct
                    cytCC.Connectivity = 8;
                    cytCC.NumObjects = size(cytCC.PixelIdxList,2);
                    cytCC.ImageSize = size(nucmask);
%                     cytstats = regionprops(cytCC,'Area');
%                     this.cellData(ti).cytMaskArea = [cytstats.Area]';
                    %to get cytMaskArea, just count the pixels in each
                    %cytoplasmic mask
                    this.cellData(ti).cytMaskArea = cellfun(@numel,cytCC.PixelIdxList)';

                    this.cellData(ti).cytLevel = zeros([numel(cytCC.PixelIdxList) numel(opts.dataChannels)]);
                    this.cellData(ti).NCratio = zeros([numel(cytCC.PixelIdxList) numel(opts.dataChannels)]);

%                     cytmask = false(size(nucmask));
%                     cytmask(cat(1,cytCC.PixelIdxList{:}))=true;
                else
                    cytCC = {};
                end

                % regionprops and indices from nuclear mask
                % follows cyt so we can shrink without introducing new var
                %-------------------------------------------

                % shrink instead of erode to preserve number of connected
                % components to not introduce mismatch between nuclear and
                % cytoplasmic mask
                if opts.nucShrinkage > 0
                    nucmask = bwmorph(nucmask,'shrink',opts.nucShrinkage);
                end

                % initialize cellData
                %---------------------------------------------------

                nucCC = bwconncomp(nucmask);
                nucstats = regionprops(nucCC, 'Area', 'Centroid','PixelIdxList',...
                                                'Orientation', 'MajorAxisLength','MinorAxisLength',...
                                                'Circularity');
                nCells = numel(nucstats);
                %align indexing between nuclear and cytoplasmic masks for
                %individual cells
                %I think cytCC should strictly go with the watershed
                %indices, so this could probably be simplified a little
                if opts.cytoplasmicLevels
                    %shuffle cytoplasmic masks to be in the same order as
                    %nuclear masks
                    cytPixelLists = cell(1,nCells);
                    cytMaskArea = NaN(nCells,1);
                    cytidxs = (1:nCells)'; nucidxs = NaN(nCells,1);
                    for cellidx = 1:nCells
                        %in theory we don't have to get the unique value
                        %over all pixels in the nucmask but it confirms
                        %that there is only one unique value in the mask; 
                        %this would result in an error otherwise
                        nucidxs(cellidx) = unique(L(nucstats(cellidx).PixelIdxList));
                        %double check that indexing for cyto pixel indices
                        %is really working the way I think
                        ind = unique(L(cytCC.PixelIdxList{cellidx}));
                        if ~isempty(ind)
                            if ind ~= cellidx
                                error('cytCC indexing is wrong')
                            end
                        end
                        
                        cytidx = find(cytidxs == nucidxs(cellidx));
                        cytPixelLists{cellidx} = cytCC.PixelIdxList{cytidx};
                        cytMaskArea(cellidx) = this.cellData(ti).cytMaskArea(cytidx);
                    end
                    cytCC.PixelIdxList = cytPixelLists;
                    this.cellData(ti).cytMaskArea = cytMaskArea;
                end
                
                
                if keepPixelIdxs
                    this.cellData(ti).nucmask = {nucstats.PixelIdxList};
                    if opts.cytoplasmicLevels
                        this.cellData(ti).cytmask = cytCC.PixelIdxList;
                    end
                end
                
                this.ncells(ti) = nCells;
                this.cellData(ti).XY = cat(1,nucstats.Centroid);
                
                % nuclear geometry
                this.cellData(ti).nucArea = cat(1,nucstats.Area);
                this.cellData(ti).nucOrientation = cat(1,nucstats.Orientation);
                this.cellData(ti).nucMajorAxis = cat(1,nucstats.MajorAxisLength);
                this.cellData(ti).nucMinorAxis = cat(1,nucstats.MinorAxisLength);
                this.cellData(ti).nucCircularity = cat(1,nucstats.Circularity);
                
                this.cellData(ti).nucLevel = zeros([nCells numel(opts.dataChannels)]);
                this.cellData(ti).nucZ = zeros([nCells 1]);

                this.cellData(ti).background = zeros([1 numel(opts.dataChannels)]);

                % if is no MIPidx, analyze the MIP
                % taking the MIP of a single z-slice costs no time so
                % this works for single images
                if ~isempty(opts.MIPidxDir)
                    MIPidx = this.loadMIPidx(opts.MIPidxDir, this.nucChannel, ti);
                else
                    MIPidx = [];
                end
                % ti == 1 so it doesnt say it a hundred times
                if ti == 1 && isempty(MIPidx)
                   warning('------------ NO MIPidx FOUND ------------');
                end

                % for background subtraction, median z-plane
                % this ifempty(MIPidx) is just so that it proceeds without
                % MIPidx if the MIPidx cannot be loaded
                if ~isempty(MIPidx)
                    zmed = median(MIPidx(nucmask));
                else
                    % if there is no MIPidx, it will use the MIP
                    % read out in its only plane
                    zmed = 1;
                end

                % read out nuclear and cytoplasmic levels
                %-----------------------------------------
                % zmed == 0 means no nuclei in mask so no cells in MIPidx
                % case
                if nCells > 0  % zmed > 0

                for cii = 1:numel(opts.dataChannels)

                    %disp(['loading channel ' num2str(opts.dataChannels(cii))]);
                    if strcmp(ext,'.oib')
                        imc = squeeze(imgdata(:,:,opts.dataChannels(cii),:,ti));
                    else
                        imc = this.loadImage(dataDir, opts.dataChannels(cii), ti);
                    end
                    %disp(['size: ' num2str(size(imc))]);

                    if size(nucmask) ~= [size(imc,1) size(imc,2)]
                        error(['nucmask size ' num2str(size(nucmask)) ' does not match image size ' num2str(size(imc))]);
                    end

                    % if no MIPidx, just use MIP
                    if isempty(MIPidx)
                        imc = max(imc,[],3);
                    end

                    if ~isempty(opts.imopenBGsub)
                        disp('imopen bg sub');
                        imc = imc - imopen(imc,strel('disk',opts.imopenBGsub));
                    end

                    % current low-tech background subtraction:
                    %-------------------------------------------------
                    % mean value of segmented empty space in the image
                    % or otherwise just min of image
                    if ~isempty(bgmask)
                        % disp('bg seg bg sub');
                        % if the background area is too small to be
                        % reliable (I'm saying < ~1% of total area), then
                        % just go with the previous value
                        % zmed = 0 means no more nuclei left
                        
                        imcZmed = imc(:,:,zmed);
                        %dont use black buffer regions around the image for
                        %background estimation
                        bgmask = bgmask & (imcZmed > 0);
                        if sum(bgmask,'all') > numel(bgmask)/100
                            if ~isempty(MIPidx)
                                %use weighted average for background estimation
                                vals = MIPidx(nucmask);
                                indices = unique(vals);
                                npoints = numel(vals);
                                bg_est = 0;
                                for jj = 1:length(indices)
                                    ind = indices(jj);
                                    imcZmed = imc(:,:,ind);
                                    weight = sum(vals == ind)/npoints;
                                    bg_est = bg_est + weight*mean(imcZmed(bgmask));
                                end
                                this.cellData(ti).background(cii) = bg_est;
                            else
                                this.cellData(ti).background(cii) = mean(imcZmed(bgmask));
                            end
                        else
                            if ti > 1
                                this.cellData(ti).background(cii) = this.cellData(ti-1).background(cii);
                            else
                                this.cellData(ti).background(cii) = min(min(imcZmed));
                            end
                        end
                    else
                        imcZmed = imc(:,:,zmed);
                        this.cellData(ti).background(cii) = min(min(imcZmed(imcZmed>0)));
                    end

                    for cellidx = 1:nCells

                        if ~isempty(MIPidx)
                            zi = median(MIPidx(nucstats(cellidx).PixelIdxList));
                        else
                            zi = 1;
                        end

                        nucPixIdx = nucstats(cellidx).PixelIdxList;
                        nucPixIdx = double(zi-1)*size(imc,1)*size(imc,2) + nucPixIdx;
                        this.cellData(ti).nucLevel(cellidx, cii) = mean(imc(nucPixIdx));
                        this.cellData(ti).nucZ(cellidx) = zi;

                        if opts.cytoplasmicLevels
                            
                            cytPixIdx = cytCC.PixelIdxList{cellidx};
                            cytPixIdx = double(zi-1)*size(imc,1)*size(imc,2) + cytPixIdx;
                            this.cellData(ti).cytLevel(cellidx, cii) = mean(imc(cytPixIdx));
                            
                            C = this.cellData(ti).cytLevel(cellidx, cii);
                            N = this.cellData(ti).nucLevel(cellidx, cii);
                            bg = this.cellData(ti).background(cii);
                            this.cellData(ti).NCratio(cellidx, cii) = (N-bg)./(C-bg);
                        end
                    end
                end
                else
                    warning(['------------ NO CELLS AT T = ' num2str(ti) '------------']);
                end
            end
            fprintf('\n');

            % make cytomask for debugInfo
            cytmask = false(size(nucmask));
            if opts.cytoplasmicLevels
                cytmask(cat(1,cytCC.PixelIdxList{:}))=true;
            end

            debugInfo = struct('bgmask', bgmask, 'nucmaskraw', nucmaskraw, 'nucmask', nucmask,...
                                'fgmask', fgmask, 'cytmask', cytmask,'segall',{seg});
        end

        function [CM_size, comp_time] = makeTracks(this, max_linking_distance,...
    max_gap_distance, max_gap_closing, imsize)

            %%%% Inputs %%%%
            % max_linking_distance      maximum distance for frame-frame cell movement
            % max_gap_distance          maximum cell movement distance when doing gap closing (distance moved over multiple frames)
            % max_gap_closing           maximum # of frames for gap-closing, merging, and splitting
            % imsize                    size of input image (get this from meta?)
            %%%% Extract data from this %%%%
            points = transpose({this.cellData.XY});
            if isfield(this.cellData,'nucArea')
                areas = transpose({this.cellData.nucArea});
            else
                areas = transpose({this.cellData.area});
            end

            intensities = transpose({this.cellData.nucLevel});
            for i=1:length(intensities)
                %get nuclear levels
                intensities{i} = intensities{i}(:,this.nucChannel + 1);
            end
            ntimepoints = this.nTime;

            %get rid of empty time points (pad with NaNs)
            emptyTimes = cellfun(@isempty, points);
            if sum(emptyTimes == 1) > 0
                for i = find(emptyTimes)
                    points{i} = [NaN NaN];
                end
            end

            %%%% Frame-frame linking %%%%
            disp('Doing frame-frame linking step')
            tic
            %Variables to hold data during linking
            target_indices = cell(size(points));
            source_indices = cell(size(points));
            %Frame-to-frame linking step
            parfor i=1:(ntimepoints - 1)
                source = [points{i}, areas{i}, intensities{i}];
                target = [points{i+1}, areas{i+1}, intensities{i+1}];
                [ target_idxs, source_idxs ] = ...
                    newhungarianlinker(source, target, imsize, max_linking_distance);
                target_indices{i} = target_idxs;
                source_indices{i} = source_idxs;
                if mod(i,round((ntimepoints-1)/20)) == 0
                    disp('.')
                end
            end
            toc

            %%%% String frame-to-frame links into tracks over the full time period %%%%
            %Simple tracks is a cell object that holds tracks with no gap closing,
            %linking or merging
            simple_tracks = cell(size(points{1},1), 1);

            %goods is used to keep track of which tracks are still active; for each
            %frame, if a track i ended (linked to -1), then goods(i) is set to zero,
            %and if a new track is started, a new array is added to simple_tracks and a
            %new value of 1 is appended to goods
            goods = ones(size(simple_tracks));
            for i=1:size(points{1},1)
                simple_tracks{i} = [i, 1];
            end

            for i=1:ntimepoints - 1
                %In each frame, only look at tracks for which goods(j) = 1; this gives
                %the indices of all such tracks in simple_tracks
                source_tracks = find(goods == 1);
                for j=1:length(source_tracks)
                    jj = source_tracks(j);
                    sp = simple_tracks{jj}(end,1);
                    tp = target_indices{i}(sp);
                    %If the point in track j is linked to nothing in the next frame,
                    %set goods(jj) = 0; this means the track ends in this frame
                    if tp == -1
                        goods(jj) = 0;
                    else
                        %Otherwise add the index of the point in the next frame to the
                        %row of simple_tracks(jj) alongside its time index
                        simple_tracks{jj} = [simple_tracks{jj};[tp, i+1]];
                    end
                end
                %Find points at time i which are linked to nothing in the previous
                %frame; these correspond to new tracks
                new_idxs = find(source_indices{i} == -1);
                for j=1:length(new_idxs)
                    jj = new_idxs(j);
                    %Add an entry to simple tracks starting with the index and time of
                    %the first point in the track
                    simple_tracks{size(simple_tracks,1)+1} = [jj, i+1];
                    %Append a value of 1 to goods to indicate that this track is active
                    goods = [goods; 1];
                end
            end

            ntracks = size(simple_tracks, 1);
            disp(strcat("Number of simple tracks = ", num2str(ntracks)))

            %%%% Gap closing, splitting, merging, creating compound tracks %%%%
            [gaps, merges, splits, ~, CM_size, comp_time] = tracklinker8(simple_tracks,...
                 max_gap_closing, max_gap_distance, points, intensities);
             %Output: nx4 arrays with linking info
            %gaps(i,:) = [end_track, start_track, gap_time, cost]
            %merges(i,:) = [end_track, mid_track, merge_time, cost]
            %splits(i,:) = [end_track, mid_track, split_time, cost]


            %%%% Close gaps in simple_tracks %%%%
            disp("Closing gaps in simple tracks")
            tCol = 2;
            idxCol = 1;
            %Start by grouping gap-closed tracks together - track_idxs should hold a
            %cell array where indices of these tracks are grouped; i.e. if track 3 is
            %closed to track 12 and the end of track 12 links to the beginning of track
            %30 and track 30 terminates without another link, then one of the entries
            %of track_idxs will be [3,12,30], and the track_ids for 12 and 30 will be
            %changed to 0
            %track_idxs{i,2} holds a nx2 array with costs and times for gap
            %closes done between the simple_tracks in track_idxs{i,1}
            track_idxs = cell(ntracks,1);
            gap_conf = cell(ntracks,1);
            track_ids = ones(ntracks,1);
            for i=1:ntracks
                track_ids(i) = 1; %unnnecessary?
                track_idxs{i,1} = i;
            end

            for i = 1:ntracks
                idx = i;
                if ismember(idx,gaps(:,1)) && track_ids(idx)
                    idxlist = idx;
                    cts = []; %costs and times
                    criterion = true;
                    while criterion == true
                        idxNext = gaps(gaps(:,1) == idx,2);
                        cost_time = gaps(gaps(:,1) == idx,3:4);
                        track_ids(idxNext) = 0;
                        idxlist = [idxlist,idxNext];
                        cts = [cts; cost_time];
                        if ismember(idxNext,gaps(:,1))
                            if idx == idxNext
                            %prevent infinite while loop when there is an error in CM construction
                                error('Incorrect track assignments')
                            end
                            idx = idxNext;
                        else
                            criterion = false;
                        end
                    end
                    track_idxs{i,1} = idxlist;
                    gap_conf{i,1} = cts;
                end
            end

            %Revise the track_idxs to get rid of duplicates using track_ids and reset
            %track_ids for merging and splitting
            track_idxs2 = cell(sum(track_ids),1);
            gap_conf2 = cell(sum(track_ids), 1);
            idx = 1;
            for i=1:size(track_idxs,1)
                if track_ids(i) == 1
                    track_idxs2{idx,1} = track_idxs{i,1}; %track indices
                    gap_conf2{idx,1} = gap_conf{i,1}; %costs for linking & times
                    idx = idx + 1;
                end
            end
            track_idxs = track_idxs2;
            gap_conf = gap_conf2;
            track_ids = ones(size(track_idxs,1),1);

            %Make gap_tracks using track_idxs; this makes groups of tracks that are
            %gap-closed into single tracks with discontinuities
            gap_tracks = cell(length(track_ids),1);
            for i=1:size(gap_tracks,1)
                %time points with no cell have nan for the cell index
                gap_tracks{i} = nan(ntimepoints,1);
                for j=1:length(track_idxs{i,1})
                    sidx = track_idxs{i,1}(j); %index of corresponding simple track
                    times = simple_tracks{sidx}(:,tCol); %times for which the track is defined
                    idxs = simple_tracks{sidx}(:,idxCol); %cell index in frame at these times
                    gap_tracks{i}(times) = idxs;
                end
            end
            %gap_conf{i,1} still has the costs and times for gap_tracks{i}

            %Go through merges and splits and replace track indices of simple_tracks
            %with corresponding track indices for gap_tracks
            splits2 = splits;
            merges2 = merges;
            for i = 1:size(splits2,1)
                for j = 1:2
                    sidx = splits2(i,j);
                    gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
                    splits2(i,j) = gidx;
                end
            end

            for i = 1:size(merges2,1)
                for j = 1:2
                    sidx = merges2(i,j);
                    gidx = find(cellfun(@(x) ismember(sidx,x), track_idxs));
                    merges2(i,j) = gidx;
                end
            end

            %%%% Handle merging events %%%%
            %Merging events are assumed to be artifacts of the detection, segmentation,
            %or tracking process and should be handled before using splitting events to
            %create compound tracks

            disp("Resolving merging and splitting events")
            %Find which tracks are merged to
            merge_tracks = unique(merges2(:,2));

            %For each track that is merged to, create a list of events that happen to
            %it: merges and splits, the tracks that are merging and splitting to/from
            %them, and the time point at which event occurs; label merges and splits
            %with 0 and 1, respectively
            %[merge_track, other track, time, merge/split (0/1)]
            events = cell(length(merge_tracks),1);
            for k = 1:length(merge_tracks)
                mevents = merges2(merges2(:,2) == merge_tracks(k),[2,1,3]);
                sevents = splits2(splits2(:,2) == merge_tracks(k),[2,1,3]);
                events{k} = [mevents, zeros(size(mevents,1),1); sevents, ones(size(sevents,1),1)];
                [~, I] = sort(events{k}(:,3));
                events{k} = events{k}(I,:);
            end

            %For each unique track that experiences at least one merging event, create
            %a cost matrix to make assignments: each of the two branches in a merging
            %event may be assigned either to a subsequent split branch or considered an
            %error
            %Set a cutoff such that if there are more than so many frames after a
            %merge before another split, it is considereed an erroneous assignment

            %For now, essentially just gap-close, assuming that if there is a split
            %soon after a merge then the merged track corresponds to the split track
            %and otherwise that the merged track terminates; this is not robust
            merge_cutoff = 10;
            [newtracks, splits4] = newmergelinker(events, gap_tracks,...
            splits2, points, areas, intensities, merge_cutoff);

            disp("Constructing final tracks")
            this.tree = build_tree(newtracks, splits4);

            nctracks = numel(getchildren(this.tree,1));
            disp(strcat("Number of compound tracks = ", num2str(nctracks)))

        end

        function makeAvgTimeTraces(this, SNRcutoff)
            % make time traces of levels in cellData
            %
            % makeTimeTraces()
            % makeTimeTraces(SNRcutoff)
            %
            % populates Position.timeTraces
            %
            % SNRcutoff :   exclude cell if N/bg or C/bg < SNRcutoff,
            %               default 2

            if ~exist('SNRcutoff','var')
                SNRcutoff = zeros([numel(this.dataChannels) 1]);
            end

            nucLevelAvg = zeros([this.nTime numel(this.dataChannels)]);
            nucLevelMed = zeros([this.nTime numel(this.dataChannels)]);

            cytLevelAvg = zeros([this.nTime numel(this.dataChannels)]);
            cytLevelMed = zeros([this.nTime numel(this.dataChannels)]);
            ratioAvg = zeros([this.nTime numel(this.dataChannels)]);

            bg = zeros([this.nTime numel(this.dataChannels)]);

            for ti = 1:numel(this.cellData)

                nucLevel = this.cellData(ti).nucLevel;

                if isfield(this.cellData(ti),'cytLevel')
                    cytLevel = this.cellData(ti).cytLevel;
                else
                    cytLevel = [];
                end

                A = this.cellData(ti).nucArea;

                for ci = 1:numel(this.dataChannels)

                    bg(ti,ci) = this.cellData(ti).background(ci);

                    % implement SNR cutoff necessary for varying expression
                    % levels to exclude outliers from dim cells
                    if ~isempty(cytLevel)
                        tooDim =    cytLevel(:,ci)/bg(ci) < SNRcutoff(ci) |...
                                    nucLevel(:,ci)/bg(ci) < SNRcutoff(ci);
                        idx = ~tooDim & ~isnan(nucLevel(:,ci)) & ~isnan(cytLevel(:,ci));
                    else
                        idx = ~isnan(nucLevel(:,ci));
                    end

                    if sum(idx)/numel(idx) < 0.2
                        warning([num2str(100*(1-sum(idx)/numel(idx)),3) '% of cells excluded based on SNR in channel ' num2str(ci) ' at time ' num2str(ti) '!']);
                    end

                    if ~isempty(idx)

                        % weight by cell area - may or may not be sensible
                        % in context
%                         W = A(idx)/mean(A(idx));
                        W = ones(size(A(idx)));
                        nucLevelAvg(ti, ci) = mean(nucLevel(idx,ci).*W);
                        nucLevelMed(ti, ci) = median(nucLevel(idx,ci).*A(idx))/mean(A(idx));

                        if ~isempty(cytLevel)
                            cytLevelAvg(ti, ci) = mean(cytLevel(idx,ci).*W);
                            cytLevelMed(ti, ci) = median(cytLevel(idx,ci).*A(idx))/mean(A(idx));
                            ratioAvg(ti,ci) = mean((nucLevel(idx,ci) - bg(ti,ci))./(cytLevel(idx,ci) - bg(ti,ci)));
                        end
                    end
                end


            end

            this.timeTraces.nucLevelAvg = nucLevelAvg;
            this.timeTraces.nucLevelMed = nucLevelMed;

            this.timeTraces.cytLevelAvg = cytLevelAvg;
            this.timeTraces.cytLevelMed = cytLevelMed;
            this.timeTraces.ratioAvg = ratioAvg;

            this.timeTraces.background = bg;
        end

        function makeTimeTraces(this, minlength)
            % this.makeTimeTraces(minlength)
            % minlength : minimal length of tracks kept (in frames)

            if ~exist('minlength','var')
                minlength = 10;
            end

            if isempty(this.tree)
                error('You must create tracks first')
            end
            
            %these are the fields (besides trackT and pointIdx) which are
            %currently saved for single cells in this.timeTraces
            fields = {'XY','NCratio','cytLevel','nucLevel','divs',...
                'nucArea', 'nucMajorAxis', 'nucMinorAxis'};
            nfields = length(fields);
            ntime = this.nTime;
            %find first frame with at least one cell
            frame = find(this.ncells > 0, 1);

            %find all the nodes of the tree with no children
            leaves = this.tree.findleaves();
            %store the total # of nodes in the tree
            nnodes = this.tree.nnodes();
            tracks = cell(size(leaves));
            %use goods to store indices of tracks with length > minlength
            goods = zeros(length(leaves),1);

            for bi = 1:length(leaves)
                %get list of node IDs for this leaf
                idx = leaves(bi);
                nodeIDs = NaN(nnodes,1);
                count = 1;
                while idx ~= 1
                    nodeIDs(count) = idx;
                    idx = this.tree.getparent(idx);
                    count = count + 1;
                end
                %put the nodes in order from root to leaf
                nodeIDs = sort(nodeIDs(~isnan(nodeIDs)));

                tracks{bi} = NaN(ntime, 1);
                for ni = 1:length(nodeIDs)
                    content = get(this.tree, nodeIDs(ni));
                    nnans = ~isnan(content);
                    tracks{bi}(nnans) = content(nnans);
                end

                %if the track is long enough, change the value of goods at
                %the index of the track to 1
                if sum(~isnan(tracks{bi})) >= minlength
                    goods(bi) = 1;
                end
            end

            %find the number of long-enough tracks, and their indices
            ntraces = sum(goods);
            goodtracks = find(goods == 1);
            %initialize cell arrays
            trackT = cell(ntraces, 1);
            pointIdx = cell(ntraces, 1);

            %initialize structs
            for fi = 1:length(fields)
                this.timeTraces.(fields{fi}) = cell(ntraces,1);
            end

            for bi = 1:length(goodtracks)
                %get track index
                tidx = goodtracks(bi);
                %get relevant time points and itialize cells for those points
                trackT{bi} = find(~isnan(tracks{tidx}));
                nt = length(trackT{bi});
                pointIdx{bi} = NaN(nt,1);

                %initialize all fields for this struct
                for fi = 1:nfields
                    width = size(this.cellData(frame).(fields{fi}),2);
                    this.timeTraces.(fields{fi}){bi} = NaN(nt,width);
                end
                %for each defined time point in the track, extract data for all
                %fields listed above
                for ti = 1:nt
                    time = trackT{bi}(ti);
                    cidx = tracks{tidx}(time);
                    pointIdx{bi}(ti) = cidx;
                    for fi = 1:nfields
                        %get value of the field in cellData at the appropriate cell
                        %at the appropriate time point and add it to the time trace
                        vals = this.cellData(time).(fields{fi})(cidx,:);
                        this.timeTraces.(fields{fi}){bi}(ti,:) = vals;
                    end
                end
            end

            this.timeTraces.trackT = trackT;
            this.timeTraces.pointIdx = pointIdx;
        end

        function cellDataOverview(this, time)

            bg = this.cellData(time).background;
            nucl = mean(this.cellData(time).nucLevel);
            area = mean(this.cellData(time).nucArea);

            disp(['number of cells: ' num2str(this.ncells(time))]);
            disp(['mean nuclear area: ' num2str(round(area)) ' pixels']);
            disp(['background: ' num2str(round(bg))]);
            disp(['mean nuclear level: ' num2str(round(nucl))]);
            if isfield(this.cellData(time),'cytLevel')
                cytl = nanmean(this.cellData(time).cytLevel);
                R = nanmean(this.cellData(time).NCratio);
                disp(['mean cytoplasmic level: ' num2str(round(cytl))]);
                disp(['N:C ' num2str(R,2)]);
                disp(['C:N ' num2str(1./R,2)]);
            end
        end

        % getter for dependent properties
        %---------------------------------

        function filename = get.filename(this)
            % filename without extension

            filename = sprintf(this.filenameFormat, this.ID-1);
        end

        function barefname = get.barefilename(this)
            % filename without extension

            barefname = this.filenameFormat(1:end-4);
        end

        % setters
        %---------------------------------

        function setID(this, ID)

            this.ID = ID;
        end
    end
end
