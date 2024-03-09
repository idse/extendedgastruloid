function [img, omeMeta] = readStack(fullfname, channels, tmax, series) 
    % read the data from a single multichannel stack
    % xyczt
    
    fullfname = char(fullfname);
            
    % load the Bio-Formats library into the MATLAB environment
    autoloadBioFormats = 1;
    status = bfCheckJavaPath(autoloadBioFormats);
    assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
        'to the static Java path or add it to the Matlab path.']);

    % initialize logging
    loci.common.DebugTools.enableLogging('INFO');
    
    % create bioformats reader for file
    disp(fullfname);
    
    tic
    r = bfGetReader(fullfname);
    if ~exist('channels','var') || isempty(channels)
        channels = 1:r.getSizeC();
    end
    if ~exist('tmax','var')
        tmax = r.getSizeT();
    end
    if exist('series','var')
        r.setSeries(series-1);
    end
    img = zeros([r.getSizeY() r.getSizeX() numel(channels) r.getSizeZ() tmax], 'uint16');
    for ti = 1:tmax
        for cii = 1:numel(channels)
            for zi = 1:r.getSizeZ()
                fprintf('.');
                img(:,:,cii,zi,ti) = bfGetPlane(r, r.getIndex(zi-1,channels(cii)-1,ti-1)+1);
            end
        end
    end
    fprintf('\n')
    omeMeta = r.getMetadataStore();
    r.close();
    toc
    %img = squeeze(img); % DON'T DO THIS, makeMIP needs xyczt even if c is
    %singleton because of the channels parameter
end
