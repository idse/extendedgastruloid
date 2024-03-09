function Ilim = makeXsection(im, fname_prefix, meta, options)
    % im: image stack
    % fname_prefix : start of filename
    % meta: metadata
    % index : y index to make cross section at
    % RGBset : channels to use for RGB image
    % tol: tolerances for channel adjustment
    % Ilim : intensity limits cell array per channel
    
    
    index = options.index;
    RGBset = options.RGBset;
    tol = options.tol;
        
    if ~isfield(options,'Ilim') || isempty(options.Ilim)
        Ilim = {};
        fillIlim = true;
    else
        Ilim = options.Ilim;
        fillIlim = false;
    end
    
    if ~isfield(options, 'direction')
        direction = 'y';
    else
        direction = options.direction;
    end
    
    if ~isfield(options,'color')
        options.color = 'RGB';
    end

    if ~isfield(options,'scalebar')
        options.scalebar = true;
        % scalebar parameters
        scalebarNpixels = uint16(50/meta.xres);
        zmarg = 40;
        xmarg = 40;
        w = 12;
    end
    
    lw = 10;
    nucChannel = meta.nucChannel+1;
    
    MIP = max(im,[],4);
    %MIP = mean(im,4);

    MIPadj = MIP; 
    for ci = 1:meta.nChannels
        if fillIlim
            Ilim{ci} = stitchedlim(MIP(:,:,ci), struct('Tol',tol(ci,:)));
        end
        MIPadj(:,:,ci) = imadjust(MIP(:,:,ci),Ilim{ci}); 
    end
    
    t = round(lw/2);
    % draw line on xy MIP to indicate cross section
    if strcmp(options.direction,'y')
        MIPadj(index-t:index+t,:,:) = max(MIPadj(:));
    elseif strcmp(options.direction,'x')
        MIPadj(:,index-t:index+t,:) = max(MIPadj(:));
    else
        error('unknown direction');
    end
    
    newnslices = round(meta.zres / meta.xres*size(im,4));
    if isempty(newnslices) || newnslices==0 
        error('newnslices is empty or zero, resolution is probably set incorrectly');
    end
    if strcmp(options.direction,'y')
        crosssize = size(im,2);
    elseif strcmp(options.direction,'x')
        crosssize = size(im,1);
    end

    scaled = zeros([newnslices crosssize],'uint16');
    for ci = 1:meta.nChannels
        if strcmp(options.direction,'y')
            imxsect = squeeze(max(im(index,:,ci,:),[],1));
        elseif strcmp(options.direction,'x')
            imxsect = squeeze(max(im(:,index,ci,:),[],2));
        end
        scaled(:,:,ci) = flipud(imresize(imxsect, [crosssize newnslices])');
        % bg subtract in cross-section
        scaled(:,:,ci) = scaled(:,:,ci) - imopen(scaled(:,:,ci), strel('disk', 50));
        scaled(:,:,ci) = imadjust(scaled(:,:,ci),Ilim{ci});
        imwrite(scaled(:,:,ci), [fname_prefix '_' meta.channelLabel{ci} '_' num2str(index([1 end])) '_' direction 'z.png']);
    end
    C1 = scaled(:,:,RGBset(1));
    if numel(RGBset)>=2
        C2 = scaled(:,:,RGBset(2));
    else
        C2 = 0*C1;
    end
    if numel(RGBset)==3
        C3 = scaled(:,:,RGBset(3));
    else
        C3 = 0*C1;
    end

    if strcmp(options.color,'RGB')
        combined = cat(3, C1, C2, C3);
    elseif strcmp(options.color, 'CMY')
        combined = cat(3, C3+C2, C1+C3, C1+C2);
    elseif strcmp(options.color, 'GM')
        combined = cat(3, C2, C1, C2);
    else
        error('invalid color option');
    end
    
    % add scalebar
    if options.scalebar
        combined(end-zmarg-w:end-zmarg, end-xmarg-scalebarNpixels:end-xmarg,:) = 2^16-1;
    end

    % save
    imwrite(combined, [fname_prefix '_' meta.channelLabel{RGBset} '_' num2str(index([1 end])) '_' options.color '_' direction 'z.png']);

    figure,
    for ci = 1:meta.nChannels

        subplot_tight(meta.nChannels + 1, 1, ci)
        imshow(scaled(:,:,ci),[])
        title(meta.channelLabel{ci})
    end

    subplot_tight(meta.nChannels + 1, 1, meta.nChannels + 1)
    xzim = combined;
    if ~isempty(nucChannel)
        xzimnuc = scaled(:,:,nucChannel);
        for i = 1:3
           xzim(:,:,i) =  xzim(:,:,i) + xzimnuc;
        end
    end
    imwrite(xzim, [fname_prefix '_' meta.channelLabel{RGBset} 'DAPI_' num2str(index([1 end])) '_' options.color '_' direction 'z.png']);
    imshow(xzim,[])

    str = [];
    for i = 1:numel(RGBset)
        str = [str meta.channelLabel{RGBset(i)} '; '];
    end
    title(str);
    saveas(gcf, [fname_prefix '_' direction 'sect_' num2str(index([1 end])) '.png']);
end