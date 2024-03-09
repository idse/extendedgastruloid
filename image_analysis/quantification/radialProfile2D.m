function result = radialProfile2D(stats, meta, condi, options)
    
    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'half')
        options.half = false;
    end

    scalefactor = round(meta.zres / meta.xres);

    % get stats
    time = 1;
    Z = stats.XY{condi,time}(:,3)*scalefactor;
    R = sqrt(sum(stats.XY{condi,time}(:,1:2).^2,2));
    nucLevel = stats.nucLevel{condi,time};

    marg = 55;
    if ~isfield(options,'Rmax')
        Rmax = ceil(max(R(:)) + marg);
    else
        Rmax = options.Rmax;
    end
    if ~isfield(options, 'Zmax')
        Zmax = ceil(max(Z(:)) + marg/2);
    else
        Zmax = options.Zmax;
    end

    % create equal volume bins to obtain density
    Nrbins = 100;%round(numel(R)/100);
    Redges = linspace(0, Rmax^2, Nrbins);
    Redges = sqrt(Redges);
    % not sure why, less Z bins needed for smooth result
    Nzbins = round(Nrbins/4);
    Zedges = linspace(0, Zmax, Nzbins); 

     % make grid of bin centers (avg of neighboring bin edges)
    reqvol = conv(Redges,[1 1]/2,'valid');
    zeqvol = conv(Zedges,[1 1]/2,'valid');
    % but make the central most bin center r=0 so there is no hole in the
    % middle
    reqvol(1) = 0;
    [Rgeqvol,Zgeqvol] = meshgrid(reqvol, zeqvol);

    % create linear grid for interpolation
    Rlin = linspace(0, Rmax, Nrbins);
    Zlin = linspace(0, Zmax, round(Nrbins*(Zmax/Rmax)));
    [Rgl,Zgl] = meshgrid(Rlin, Zlin);

    nc = meta.nChannels;
    nuc_profile_tmp = NaN([Nzbins-1 Nrbins-1 nc]);
    nuc_profile_r_tmp = NaN([Nrbins-1 nc]);
    nuc_profile_z_tmp = NaN([Nzbins-1 nc]);
%     nuc_profile_std_tmp = zeros([N-1 nc]);
%     
%     cyt_profile_tmp = zeros([N-1 nc]);
%     cyt_profile_std_tmp = zeros([N-1 nc]);
%     ratio_profile_tmp = zeros([N-1 nc]);
%     ratio_profile_std_tmp = zeros([N-1 nc]);

    for i = 1:Nzbins-1

        zidx = (Z >= Zedges(i) & Z < Zedges(i+1));
        nuc_profile_z_tmp(i, :) = nanmean(nucLevel(zidx,:));
       
        for j = 1:Nrbins-1
            
            ridx = (R >= Redges(j) & R < Redges(j+1));
            combidx = zidx & ridx;
            
            if i == 1
                nuc_profile_r_tmp(j, :) = nanmean(nucLevel(ridx,:));
            end
            
            if(sum(combidx)) > 5
                nuc_profile_tmp(i,j,:) = nanmean(nucLevel(combidx,:));
            end
        end
    end
                
    % interpolate on linear grid
    nuc_profile = zeros([size(Rgl) nc]);
    nuc_profile_r = zeros([size(Rgl,2) nc]);
    nuc_profile_z = zeros([size(Rgl,1) nc]);
    
    for ci = 1:4
        nucprof_ci = interp2(Rgeqvol,Zgeqvol,nuc_profile_tmp(:,:,ci),Rgl,Zgl);
        nucprof_ci(isnan(nucprof_ci)) = min(nucprof_ci(:));
        nucprof_ci = nucprof_ci - min(nucprof_ci(:));
        nuc_profile(:,:,ci) = nucprof_ci;
        
        F=griddedInterpolant(Rgeqvol(1,:),nuc_profile_r_tmp(:,ci),'linear','nearest');
        nuc_profile_r(:,ci) = F(Rgl(1,:));
        F=griddedInterpolant(Zgeqvol(:,1),nuc_profile_z_tmp(:,ci),'linear','nearest');
        nuc_profile_z(:,ci) = F(Zgl(:,1));
    end
    
    % normalize DAPI intensity for intensity loss in z
    % only DAPI suffers from that significantly and is only needed here to
    % define to colony borders anyway
    DAPI_profile = nuc_profile(:,:,1);
    DAPI_profile = DAPI_profile./nuc_profile_z(:,1);
    DAPI_profile(isnan(DAPI_profile)) = 0;
    nuc_profile(:,:,1) = DAPI_profile;
    
    % mirror 
    % this basically provides a mirror boundary condition just at R=0 for
    % smoothing
    nuc_profile = cat(2, fliplr(nuc_profile), nuc_profile);
        
    % smoothing 
    sig = 7; % smoothing in micron
    histres = (Rlin(2)-Rlin(1))*meta.xres; %um/pixel in histogram
    sig = sig/histres;
    k = fspecial('gaussian', round(5*[sig sig]), sig);
    %k = fspecial('disk', round(sig));
    padding = 'symmetric';
    for ci = 1:4
        nuc_profile(:,:,ci) = imfilter(nuc_profile(:,:,ci), k, padding);
    end
    
    % cut back in half after smoothing if needed
    if options.half == true
        s = size(nuc_profile,2)/2;
        nuc_profile = nuc_profile(:,s:end,:);
    end

    % scale to original proportions for scatter overlays
    if options.half == true
        newsize = [Zmax Rmax];
    else
        newsize = [Zmax 2*Rmax];
    end
    nuc_profile_resized = zeros([newsize size(nuc_profile,3)]);
    for ci = 1:4
        nuc_profile_resized(:,:,ci) = imresize(nuc_profile(:,:,ci), newsize);
    end

    % collect results
    result = struct('nuc_profile', nuc_profile_resized,...
                    'R',R,'Z',Z,...
                    'Rmax',Rmax,...
                    'Zmax',Zmax);
end
