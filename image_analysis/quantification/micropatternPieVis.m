function Ilim = micropatternPieVis(dataDir, position, options)

    if ~isfield(options,'tol')
        tol = cat(2, [1 1 1 1]'*0.01, [1 1 1 1]'*0.99);
    else
        tol = options.tol;
    end
    if ~isfield(options,'channels')
        channels = 2:4;
    else
        channels = options.channels;
    end
    if ~isfield(options,'pieOrder')
        order = 1:numel(channels);
    else
        order = options.pieOrder;
        if numel(order) > numel(channels)
            error('numel(order) cannot be greater than numel(channels)');
        end
    end
    if ~isfield(options,'scalebar')
        options.scalebar = true;
    end
    if ~isfield(options,'segOverlay')
        options.segOverlay = false;
    end
    if ~isfield(options,'color')
        options.color = 'RGB';
    end
    % margin outside the disk shaped mask in pixels
    if isfield(options,'outsideMargin') 
        outsidemargin = options.outsideMargin;
    else
        outsidemargin = 150; 
    end
    % margin in disk shaped mask of colony in micron
    if ~isfield(options,'radialMargin') 
        options.radialMargin = 25;
    end
    
    margin = 300;
    PI = atan(1)*4;

    % self stitched files start with stitched
    if strcmp(position.filenameFormat(1:8), 'stitched')

        s = strsplit(position.filename,{'_'});
        
        if ~isempty(dir(fullfile(dataDir,'filtered',['stitched_filtered_MIP_' s{2} '*'])))
            disp('using filtered MIP')
            prefix = [s{1} '_filtered_MIP_' s{2}];
            subDir = 'filtered';
        else
            prefix = [s{1} '_MIP_' s{2}];
            subDir = 'MIP';
        end
        fnameFormat = [prefix '_w%.4d_t0000.tif'];
        
    % previous case - DOCUMENT WHAT THIS IS EXPECTING
    else
        s = strsplit(position.filename,{'_FusionStitcher','_Stitched','.ims','.tif'});
        prefix = [s{:}];
        subDir = [prefix '_zslices'];
        fnameFormat = [prefix '_MIP_w%.4d.tif'];
    end
    
    imgs = {};
    Ilim = {};
    for cii = 1:4
        
        fullfname = fullfile(dataDir, subDir, sprintf(fnameFormat, cii-1));
        try
            img = imread(fullfname);
        catch
            img = imread([fullfname(1:end-4) '.jpg']);
        end
        
        if ~isfield(options,'Ilim')
            Ilim{cii} = stretchlim(img,tol(cii,:));
        else
            Ilim{cii} = options.Ilim{cii};
        end
        imgs{cii} = imadjust(img,Ilim{cii});

        %add margin
        X = zeros(size(img)+2*margin);
        X(margin+1:margin+size(img,1),margin+1:margin+size(img,2)) = imgs{cii};
        imgs{cii} = X;
    end

    xres = position.radiusMicron / position.radiusPixel;
    center = position.center + margin;
    Rmax = uint16(position.radiusPixel + options.radialMargin/xres); % micron margin
    [X,Y] = meshgrid(1:size(imgs{1},2),1:size(imgs{1},1));
    R = sqrt((X - center(1)).^2 + (Y - center(2)).^2);
    F = atan2((X - center(1)),-(Y - center(2)));
    disk = R > Rmax;

    imgs_tmp = imgs;
    for cii = 1:4
        imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
        imgs_tmp{cii}(disk) = 1; % change 1 to 0 to get black background
    end
    
    % centered crop indices
    CM = uint16(center);
    Rcrop = Rmax + round(outsidemargin/2);
    yrange = CM(2)-Rcrop:CM(2)+Rcrop;
    xrange = CM(1)-Rcrop:CM(1)+Rcrop;
    
    % 50um scalebar on DAPI image
    Npixels = uint16(50/xres);
    marg = 50;
    w = 8;
    DAPIim = mat2gray(imgs_tmp{1}(yrange,xrange));
    if options.scalebar
        DAPIim(end-marg:end-marg+w, marg:marg+Npixels) = 0;
    end
    imwrite(DAPIim, fullfile(dataDir, sprintf([prefix '_DAPI.png'])));
        
    %-----------------------
    if options.segOverlay
        fname = sprintf([position.filename(1:end-4) '_SegOverlay.png']);
        segOverlay_in = imread(fullfile(dataDir, fname));
        X = zeros([size(img)+2*margin 3],'uint8');
        X(margin+1:margin+size(img,1),margin+1:margin+size(img,2),:) = segOverlay_in;
        segOverlay_out = X;
        for i = 1:3
            im = segOverlay_out(:,:,i);
            im(disk)=255;
            segOverlay_out(:,:,i) = im;
        end
        imwrite(segOverlay_out(yrange,xrange,:), fullfile(dataDir, sprintf([prefix '_SegOverlay_centered.png'])));
    end
    %-----------------------
    
    imgs_tmp = imgs_tmp(channels);
    if numel(channels) <= 2
        imgs_tmp{3} = disk;
    end
    if numel(channels) == 1
        imgs_tmp{2} = disk;
    end
    if strcmp(options.color,'CMY')
        overlay = cat(3, imgs_tmp{2} + imgs_tmp{3}, imgs_tmp{1} + imgs_tmp{3}, imgs_tmp{1} + imgs_tmp{2});
    elseif strcmp(options.color,'RGB')
        overlay = cat(3,imgs_tmp{:});
    elseif strcmp(options.color,'GM')
        overlay = cat(3, imgs_tmp{2}, imgs_tmp{1}, imgs_tmp{2});
    else
        error('invalid color option');
    end
        
    overlay = overlay(yrange,xrange,:);
    if options.scalebar
        overlay(end-marg:end-marg+w, marg:marg+Npixels,:) = 0;
    end
    imwrite(overlay, fullfile(dataDir, sprintf([prefix '_combined' options.color '.png'])));
    imwrite(0.8*overlay + 0.8*DAPIim, fullfile(dataDir, sprintf([prefix '_combined' options.color 'DAPI.png'])));
    
    for ci = 1:numel(channels)
        imwrite(imgs_tmp{ci}(yrange,xrange), fullfile(dataDir, sprintf([prefix '_w%.4d.png'],ci)));
    end
    
    %-----------------------
    mask = {}; 
    N = numel(channels);
    for i = 1:N
        mask{i} = F < -PI + (i-1)*2*PI/N | F > -PI + i*2*PI/N;
        mask{i} = imdilate(mask{i},strel('disk',10));
    end
    lines = mask{1};
    for i = 2:N
        lines = lines & mask{i};
    end

    imgs_tmp = imgs_tmp(order);
    for cii = 1:numel(channels)
        imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
        imgs_tmp{cii}(mask{cii}) = 0;
    end
    combined = sum(cat(3,imgs_tmp{:}),3);
    combined(disk) = 1;
    combined(lines) = 1;
    combined = combined(yrange,xrange,:);
    %imshow(combined,[])
    if options.scalebar
        combined(end-marg:end-marg+w, marg:marg+Npixels) = 0;
    end
    imwrite(combined, fullfile(dataDir, sprintf([prefix '_combined.png'])));

    % % color version
    % 
    % imgs_tmp = imgs;
    % for cii = 1:numel(channels)
    %     imgs_tmp{cii} = mat2gray(imgs_tmp{cii});
    %     imgs_tmp{cii}(mask{cii}) = 0;
    %     imgs_tmp{cii}(disk) = 1;
    %     imgs_tmp{cii}(lines) = 1;
    % end
    % combined = cat(3,imgs_tmp{:});
    % imshow(combined,[])
    
end