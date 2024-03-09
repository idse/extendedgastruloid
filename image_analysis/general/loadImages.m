function img = loadImages(fnameformat, dataDir, ID, channels, time)
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

if ~exist('time','var')
    time = 1;
end

[~,~,ext] = fileparts(fnameformat);

if strcmp(ext,'.tif') || strcmp(ext,'.btf')

    s = strsplit(fnameformat,'_F[0-9]','DelimiterType','RegularExpression');
    if numel(s) > 1
        fnameFormat = [s{1} '_MIP_p%.4d_w%.4d.tif'];
    else
        fnameFormat = fnameformat;
    end

    nchan = numel(channels);
    img = cell(1,nchan);
    for cii = 1:nchan

        if contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
            fname = sprintf(ID - 1, channels(cii), time-1);
            iminf = imfinfo(fullfile(dataDir,fname));
            im = zeros([iminf(1).Height iminf(1).Width numel(iminf)],'uint16');
            for zii = 1:numel(iminf)
                im(:,:,zii) = imread(fullfile(dataDir, fname), zii);
            end
            img{cii} = im;

        elseif ~contains(fnameFormat,'_t') && (contains(fnameFormat,'_c') || contains(fnameFormat,'_w'))
            % this doesn't handle z and t in the same tiff
            fname = sprintf(fnameFormat,ID - 1, channels(cii));
            img{cii} = imread(fullfile(dataDir, fname), time);
        else
            % this doesn't handle z and t in the same tiff
            fname = sprintf(fnameFormat,ID - 1);
            img{cii} = imread(fullfile(dataDir, fname), time);
        end
    end
    img = cat(3,img{:});

% bioformats
elseif strcmp(ext,'.vsi') || strcmp(ext, '.oif') ...
        || strcmp(ext, '.oib') || strcmp(ext, '.ims') || strcmp(ext, '.nd2')

    fname = fullfile(dataDir, sprintf(fnameformat,ID-1));

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



