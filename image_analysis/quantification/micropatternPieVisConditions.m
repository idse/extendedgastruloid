function micropatternPieVisConditions(dataDir, positions, options, meta)
    % options = struct();
    % options.channels = [2 3 4];
    % tolerances in original order of channels
    % options.tol = [0.01 0.99; 0.1 0.99; 0.01 0.99; 0.7 0.995];
    % options.positionIdx = [1 5 9]; 
    % options.invisiblePos = 23; for contrast settings 
    % options.scalebar = false;
    % 
    % micropatternPieVisConditions(dataDir, positions, options, meta)
    % 
    % dataDir can be cell array if different dirs for different pos
    
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
    if isfield(options,'invisiblePos')
        M = numel(options.invisiblePos);
        if M > 0
            invisiblePositions = positions(options.invisiblePos);
            disp('using invisible positions for contrast');
        end
    else
        M = 0;
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
    if ~isfield(options,'scalebar')
        options.scalebar = false;
    end
    
    % determine condition index for each position
    conditionIdx = [];
    for pi = options.positionIdx
        condi = find(pi >= meta.conditionStartPos,1,'last');
        conditionIdx = [conditionIdx condi];
    end
    conditionIdx = unique(conditionIdx);
    positions = positions(options.positionIdx);
    
    PI = atan(1)*4;
    margin = 300;
    N = numel(positions);
    
    % read data and find contrast limits
    imgs = {};
    Ilim = {};
    for pi = 1:N
    
%         %s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
%         s = strsplit(positions(pi).filename,{'_FusionStitcher','_Stitched','.ims','.tif'});
%         prefix = [s{:}];
%         subDir = [prefix '_zslices'];

        disp(['position ' num2str(pi)]);

        % self stitched files start with stitched
        if strcmp(positions(pi).filenameFormat(1:8), 'stitched')

            s = strsplit(positions(pi).filename,{'_'});

            % check if filtered version is present
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
            s = strsplit(positions(pi).filename,{'_FusionStitcher','_Stitched','.ims','.tif'});
            prefix = [s{:}];
            subDir = [prefix '_zslices'];
            fnameFormat = [prefix '_MIP_w%.4d.tif'];
        end

        if iscell(dataDir)
            d = dataDir{pi};
        else
            d = dataDir;
        end
        
        for cii = 1:4
            
            fname = fullfile(d, subDir, sprintf(fnameFormat, cii-1));
            try
                img = imread(fname);
            catch
                img = max(positions(pi).loadImage(d,cii-1,1),[],3);
            end
            imgs{pi, cii} = img;
            lowhigh = stretchlim(img,tol(cii,:));
            if pi == 1
                Ilim{cii} = lowhigh;
            else
                Ilim{cii} = [min(Ilim{cii}(1), lowhigh(1)) max(Ilim{cii}(2), lowhigh(2))];
            end
        end
    end
    
    % read additional positions for contrast
    for pi = 1:M
    
        s = strsplit(invisiblePositions(pi).filename,{'_FusionStitcher','.ims'});
        prefix = [s{:}];
        subDir = [prefix '_zslices'];
        if iscell(dataDir)
            d = dataDir{pi};
        else
            d = dataDir;
        end
        
        for cii = 1:4
            fname = fullfile(d, subDir, sprintf([prefix '_MIP_w%.4d.tif'], cii-1));
            img = imread(fname);
            lowhigh = stretchlim(img,tol(cii,:));
            if pi == 1
                Ilim{cii} = lowhigh; % uncomment this to only use invisible pos
            end
            Ilim{cii} = [min(Ilim{cii}(1), lowhigh(1)) max(Ilim{cii}(2), lowhigh(2))];
        end
    end
    
    % adjust contrast, center and crop
    for pi = 1:N
        for cii = 1:4
            imgs{pi, cii} = imadjust(imgs{pi, cii}, Ilim{cii});
            imsize = size(imgs{pi, cii});
            %add margin
            X = zeros(size(img)+2*margin);
            X(margin+1:margin+imsize(1),margin+1:margin+imsize(2)) = imgs{pi, cii};
            imgs{pi, cii} = X;
        end

        xres = positions(pi).radiusMicron / positions(pi).radiusPixel;
        center = positions(pi).center + margin;
        Rmax = uint16(positions(pi).radiusPixel + options.radialMargin/xres); % micron margin
        
        [X,Y] = meshgrid(1:size(imgs{pi, 1},2),1:size(imgs{pi, 1},1));
        R = sqrt((X - center(1)).^2 + (Y - center(2)).^2);
        F = atan2((X - center(1)),-(Y - center(2)));
        disk = R > Rmax;

        % centered crop indices
        CM = center;
        Rcrop = Rmax + round(outsidemargin/2);
        yrange = CM(2)-Rcrop:CM(2)+Rcrop;
        xrange = CM(1)-Rcrop:CM(1)+Rcrop;

        for cii = 1:4
            imgs{pi, cii} = mat2gray(imgs{pi, cii});
            imgs{pi, cii} = imgs{pi, cii}(yrange, xrange);
        end
        F = F(yrange, xrange);
        disk = disk(yrange, xrange);
    end

    %----------------------- PIE SLICE MASKS
    mask = {}; 
    for i = 1:N
        mask{i} = F < -PI + (i-1)*2*PI/N | F > -PI + i*2*PI/N;
        mask{i} = imdilate(mask{i},strel('disk',10));
    end
    lines = mask{1};
    for i = 2:N
        lines = lines & mask{i};
    end
    %-----------------------
    
    pieslices = imgs;
    for pi = 1:N
        for cii = 1:4
            pieslices{pi, cii}(mask{pi}) = 0;
            pieslices{pi, cii}(disk) = 1;
            pieslices{pi, cii}(lines) = 1;
        end
    end
    
    condstr = strcat(meta.conditions(conditionIdx),{'_'});
    condstr = [condstr{:}];
    
    combined = {};
    for cii = 1:4
        combined{cii} = sum(cat(3,pieslices{:,cii}),3);
        if iscell(dataDir)
            d = dataDir{1};
        else
            d = dataDir;
        end
        % save first position to remember where it starts
        %imwrite(pieslices{1,cii}, fullfile(d, sprintf(['poscombined_' meta.conditions{conditionIdx(1)} '_' num2str(options.positionIdx(1)) '_w%.4d.png'],cii)));
        % save combination
        imwrite(combined{cii}, fullfile(d, sprintf(['poscombined_' condstr '_' num2str(options.positionIdx) '_w%.4d.png'],cii)));
    end
    %combinedRGB = cat(3,combined{channels});
    
    % scalebar parameters
    scalebarNpixels = uint16(50/xres);
    marg = 50;
    w = 8;
    
    imgs_tmp = {};
    for i=1:3
        if i <= numel(channels)
            imgs_tmp{i} = combined{channels(i)};
        else
            imgs_tmp{i} = lines | disk;
        end
    end
    DAPIim = combined{1};
    if strcmp(options.color,'CMY')
        overlay = cat(3, imgs_tmp{2} + imgs_tmp{3}, imgs_tmp{1} + imgs_tmp{3}, imgs_tmp{1} + imgs_tmp{2});
        if options.scalebar
            overlay(end-marg:end-marg+w, marg:marg+scalebarNpixels,:) = 0;
        end
        imwrite(overlay, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_CMY.png']));
        imwrite(0.8*overlay + 0.8*DAPIim, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_DAPI_CMY.png']));
        
    elseif strcmp(options.color,'RGB')
        overlay = cat(3,imgs_tmp{:});
        if options.scalebar
            overlay(end-marg:end-marg+w, marg:marg+scalebarNpixels,:) = 0;
        end
        imwrite(overlay, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_RGB.png']));
        imwrite(0.8*overlay + 0.8*DAPIim, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_DAPI_RGB.png']));
        
    elseif strcmp(options.color,'GM')
        overlay = cat(3, imgs_tmp{2}, imgs_tmp{1}, imgs_tmp{2});
        if options.scalebar
            overlay(end-marg:end-marg+w, marg:marg+scalebarNpixels,:) = 0;
        end
        imwrite(overlay, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_MG.png']));
        imwrite(0.8*overlay + 0.8*DAPIim, fullfile(d, ['poscombined_' condstr '_' num2str(options.positionIdx) '_DAPI_MG.png']));
    else
        error('invalid color option');
    end
end