function cleanupStack(dataDir, filename, zcorrection, zcorrfactor)
    % background subtraction, intensity correction and median filter per z-slice

    % piece for filenames from micropatternPieVis

    % self stitched files start with stitched
    if strcmp(filename(1:8), 'stitched')
        s = strsplit(filename,{'_'});
        prefix = 'stitched';
        postfix = [s{2} '_w%.4d_t0000.tif'];
        fnameFormat = [prefix '_' postfix];
        cleanfnameFormat = [prefix '_filtered_' postfix];
        cleanMIPfnameFormat = [prefix '_filtered_MIP_' postfix];
        
    % previous case - DOCUMENT WHAT THIS IS EXPECTING
    else
        s = strsplit(filename,{'_FusionStitcher','_Stitched','.ims','.tif'});
        prefix = [s{:}];
        subDir = [prefix '_zslices'];
        fnameFormat = [prefix '_MIP_w%.4d.tif'];
    end
    
    % make directory for output
    if ~exist(fullfile(dataDir,'filtered'),'dir')
        mkdir(fullfile(dataDir,'filtered'));
    end
    
    if ~exist('zcorrection','var')
        zcorrection = false;
    end
    if ~exist('zcorrfactor','var')
        zcorrfactor = [1 1 1 1]*1.05;
    end
    
    imgs = {};
    Ilim = {};
    
    for channelIdx = 0:3 

        fname = sprintf(fnameFormat, channelIdx);
        fullfname = fullfile(dataDir,fname);
        
        if exist(fullfname, 'file')
            
        img = readStack(fullfname);
        img=squeeze(img);

        medstack = 0*img;
        for zidx = 1:size(img,3)
            
            % bg subtraction
            if zidx == 1
                mode = 'overwrite';
            else
                mode = 'append';
            end
            slice = img(:,:,zidx); 
            slicebg = imopen(slice,strel('disk',50));
            slicebgsub = slice-slicebg;
            
            % correct intensity
            if zcorrection
                slicebgsub = slicebgsub*(zcorrfactor(channelIdx+1)^(zidx-1));
            end
            
            %openstack(:,:,zidx) = imopen(slicebgsub,strel('disk',2));
            medstack(:,:,zidx) = medfilt2(slicebgsub,[5 5]);
            tmp = medstack(:,:,zidx);
            
            cleanfname = sprintf(cleanfnameFormat, channelIdx);
            imwrite(tmp,fullfile(dataDir,'filtered',cleanfname),'WriteMode',mode);
        end

        medMIP = max(medstack,[],3);
        cleanMIPfname = sprintf(cleanMIPfnameFormat, channelIdx);
        imwrite(medMIP,fullfile(dataDir,'filtered',cleanMIPfname));
        
        end
    end
end
