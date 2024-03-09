function result = radialPositive2D(positions, meta, channelThresholds, options)
    
    if ~exist('options','var')
        options = struct();
    end
    if ~isfield(options,'half')
        options.half = false;
    end
    
    if ~isfield(options,'channels')
        channels = [2 3 4]; 
    else
        channels = options.channels;
    end
    
    scalefactor = round(meta.zres / meta.xres);

    % combine data from different positions (if any)
    Z = [];
    R = [];
    nucLevel = [];
    for pi = 1:numel(positions)

        Zp = positions(pi).cellData.Z*scalefactor;
        Rp = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));
        nucLevelp = positions(pi).cellData.nucLevel - positions(pi).cellData.background;

        Z = cat(1, Z, Zp);
        R = cat(1, R, Rp);
        nucLevel = cat(1, nucLevel, nucLevelp);
    end

    % index for positive cells
    positive = {};
    for ci = 1:meta.nChannels
        positive{ci} = nucLevel(:,ci) > channelThresholds(ci);
    end

    marg = 55;
    if ~isfield(options,'Rmax')
        Rmax = max(R(:)) + marg;
    else
        Rmax = options.Rmax;
    end
    if ~isfield(options, 'Zmax')
        Zmax = max(Z(:)) + round(marg/2);
    else
        Zmax = options.Zmax;
    end

    % create equal volume bins to obtain density
    Nbins = round(numel(R)/100);
    Redges = linspace(0, Rmax^2, Nbins);
    Redges = sqrt(Redges);
    % not sure why, less Z bins needed for smooth result
    Zedges = linspace(0, Zmax, round(Nbins/2)); 
    totaldensity = histcounts2(R, Z, Redges, Zedges)';
    % SHOULD I PUT UNITS ON THIS? 
    
    % make grid of bin centers (avg of neighboring bin edges)
    reqvol = conv(Redges,[1 1]/2,'valid');
    zeqvol = conv(Zedges,[1 1]/2,'valid');
    % but make the central most bin center r=0 so there is no hole in the
    % middle
    reqvol(1) = 0;
    [Rgeqvol,Zgeqvol] = meshgrid(reqvol, zeqvol);
    
    % interpolate on linear grid 
    Rlin = linspace(0, Rmax, Nbins);
    Zlin = linspace(0, Zmax, round(Nbins*Zmax/Rmax));
    [Rgl,Zgl] = meshgrid(Rlin, Zlin);
    
    totaldensity = interp2(Rgeqvol,Zgeqvol,totaldensity,Rgl,Zgl);
    totaldensity(isnan(totaldensity)) = 0;
    
    %  make bins with a fixed number of cells per bin along z or r
    % (fixed number of cell/bin in marginals)
    % this prevents artifacts in low density regions like to colony edge
    ptsperbinindividualcolonies = 100;
    [Rs,Ir] = sort(R);
    %Nbins = round(numel(Ir)/ptsperbinindividualcolonies);
    Redges = unique(Rs(1:ptsperbinindividualcolonies:end));
    
    % unique is added because we ran into a case where many cells had
    % identical z values making multiple z-edges identical which causes
    % problems
    zptsperbinindividualcolonies = round(ptsperbinindividualcolonies*Rmax/Zmax);
    [Zs,Iz] = sort(Z);
    Zedges = unique(Zs(1:zptsperbinindividualcolonies:end));
    
    % combined histogram cells
    % unlike totaldensity,
    % this DOES NOT REPRESENT DENSITY since bin sizes are not fixed and in
    % fact designed to get similar counts everywhere
    % dividing by bin width would give density but not corrected for
    % geometry (less cells at smaller r)
    % this doesn't matter since we use this to normalize cells that are
    % positive or negative for markers
    totalhist = histcounts2(R, Z, Redges, Zedges)';
    
%     % alternative approach to constucting density
%     dR = Redges(2:end) - Redges(1:end-1);
%     dZ = Zedges(2:end) - Zedges(1:end-1);
%     totaldensity = totalhist./(dZ*sqrt(dR'));
    
    % any positive hist
    usedchannels = channels(~isnan(channels));
    anypositive = positive{usedchannels(1)};
    for i = 2:numel(usedchannels)
        anypositive = anypositive | positive{usedchannels(i)};
    end
    anyposhist = histcounts2(R(anypositive), Z(anypositive), Redges, Zedges)';
    anyposhist_norm = anyposhist./totalhist;
    anyposhist_norm(isnan(anyposhist_norm))=0;
    
    % negative hist
    neghist = totalhist-anyposhist;
    neghistR_norm = sum(neghist,1)./sum(totalhist,1); 
    neghistZ_norm = sum(neghist,2)./sum(totalhist,2);
    neghist_norm = neghist./totalhist;
    
    % positive for something 
    poshist = {};
    poshist_norm = {};
    poshistR_norm = {};
    poshistZ_norm = {};
    for ci = 1:3
        poshist{ci} = histcounts2(R(positive{ci+1}), Z(positive{ci+1}), Redges, Zedges)';

        % marginals
        poshistR_norm{ci} = sum(poshist{ci},1)./sum(totalhist,1);
        poshistZ_norm{ci} = sum(poshist{ci},2)./sum(totalhist,2);
        
        % normalize
        poshist_norm{ci} = poshist{ci}./totalhist;
        poshist_norm{ci}(isnan(poshist_norm{ci})) = 0;
    end

    % surf(Rg, Zg, densityhistcomb);
    % colormap gray
    % axis equal
    % caxis([0 1]);
    % axis off
    % view(2);
    
    % make grid of bin centers (avg of neighboring bin edges)
    r = conv(Redges,[1 1]/2,'valid');
    z = conv(Zedges,[1 1]/2,'valid');
    % but make the central most bin center r=0 so there is no hole in the
    % middle
    r(1) = 0;
    [Rg,Zg] = meshgrid(r, z);
    
    %----------------------------------------------------------------------
    % interpolate everything on linear grid
    %
    % normalized quantities will have NaNs where total was zero - exclude
    % these from the interpolation by using scatteredInterpolant of non-NaN
    % values
    %
    % nearest extrapolation to make sure the totals add up to 1 even on the
    % edge (otherwise we are smoothing to zero for all)
    % (they add up to 1 elsewhere so this way also in the extrapolation)
    %----------------------------------------------------------------------
    
    totalhist = interp2(Rg,Zg,totalhist,Rgl,Zgl);
    totalhist(isnan(totalhist)) = 0;
    
    F = scatteredInterpolant(Rg(:),Zg(:),neghist(:),'linear','nearest');
    neghist = F(Rgl,Zgl);

    good = ~isnan(neghist_norm(:));
    F = scatteredInterpolant(Rg(good),Zg(good),neghist_norm(good),'linear','nearest');
    neghist_norm = F(Rgl,Zgl);
  
    good = ~isnan(anyposhist_norm(:));
    F = scatteredInterpolant(Rg(good),Zg(good),anyposhist_norm(good),'linear','nearest');
    anyposhist_norm = F(Rgl,Zgl);
    
    F=griddedInterpolant(Rg(1,:),neghistR_norm,'linear','nearest');
    neghistR_norm = F(Rgl(1,:));
    F=griddedInterpolant(Zg(:,1),neghistZ_norm,'linear','nearest');
    neghistZ_norm = F(Zgl(:,1));

    for ci = 1:3
        F = scatteredInterpolant(Rg(:),Zg(:),poshist{ci}(:),'linear','nearest');
        poshist{ci} = F(Rgl,Zgl);

        good = ~isnan(poshist_norm{ci});
        F = scatteredInterpolant(Rg(good),Zg(good),poshist_norm{ci}(good),'linear','nearest');
        poshist_norm{ci} = F(Rgl,Zgl);

        F=griddedInterpolant(Rg(1,:),poshistR_norm{ci},'linear','nearest');
        poshistR_norm{ci} = F(Rgl(1,:));
        
        F=griddedInterpolant(Zg(:,1),poshistZ_norm{ci},'linear','nearest');
        poshistZ_norm{ci} = F(Zgl(:,1));
    end

    % mirror 
    % this basically provides a mirror boundary condition just at R=0 for
    % smoothing
    totaldensity = cat(2, fliplr(totaldensity), totaldensity);
    totalhist = cat(2, fliplr(totalhist), totalhist);
    neghist_norm = cat(2, fliplr(neghist_norm), neghist_norm);
    anyposhist_norm = cat(2, fliplr(anyposhist_norm), anyposhist_norm);
    neghistR_norm = cat(2, fliplr(neghistR_norm), neghistR_norm);
    neghist = cat(2, fliplr(neghist), neghist);
    for ci = 1:3
        poshist_norm{ci} = cat(2, fliplr(poshist_norm{ci}), poshist_norm{ci});
        poshistR_norm{ci} = cat(2, fliplr(poshistR_norm{ci}), poshistR_norm{ci});
        poshist{ci} = cat(2, fliplr(poshist{ci}), poshist{ci});
    end

    % smoothing 
    sig = 7; % smoothing in micron
    histres = (Rlin(2)-Rlin(1))*meta.xres; %um/pixel in histogram
    sig = sig/histres;
    k = fspecial('gaussian', round(5*[sig sig]), sig);
    %k = fspecial('disk', round(sig));
    padding = 'symmetric'; % for z-profiles to not have artefects from smoothing at z=0;

    totaldensity = imfilter(totaldensity, k, padding);
    totalhist = imfilter(totalhist, k, padding);
    neghist = imfilter(neghist, k, padding);
    neghist_norm = imfilter(neghist_norm, k, padding);
    anyposhist_norm = imfilter(anyposhist_norm, k, padding);
    neghistR_norm = imfilter(neghistR_norm, sum(k,1), padding);
    neghistZ_norm = imfilter(neghistZ_norm, sum(k,2), padding);
    for ci = 1:3
        poshist{ci} = imfilter(poshist{ci}, k, padding);
        poshist_norm{ci} = imfilter(poshist_norm{ci}, k, padding);
        poshistR_norm{ci} = imfilter(poshistR_norm{ci}, sum(k,1), padding);
        poshistZ_norm{ci} = imfilter(poshistZ_norm{ci}, sum(k,2), padding);
    end
    
    % cut back in half after smoothing if needed
    if options.half == true
        s = size(neghist,2)/2+1;
        neghist = neghist(:,s:end);
        neghist_norm = neghist_norm(:,s:end);
        anyposhist_norm = anyposhist_norm(:,s:end);
        neghistR_norm = neghistR_norm(s:end);
        totaldensity = totaldensity(:,s:end);
        totalhist = totalhist(:,s:end);
        for ci=1:3
            poshist{ci} = poshist{ci}(:,s:end);
            poshist_norm{ci} = poshist_norm{ci}(:,s:end);
            poshistR_norm{ci} = poshistR_norm{ci}(s:end);
        end
    end
    
    % density outline: half-median density cutoff
    %totaldensity = anyposhist_norm + neghist_norm;
    %totaldensity = totalhist;
    %densitycutoff = median(totaldensity(:));
    densitycutoff = max(totaldensity(:))*0.1;
    mask = totaldensity > densitycutoff;

    % scale to original proportions for scatter overlays
    if options.half == true
        newsize = [Zmax Rmax];
    else
        newsize = [Zmax 2*Rmax];
    end
    density_scaled = imresize(totaldensity, newsize);
    mask_scaled = density_scaled > densitycutoff; %mask_scaled = imresize(mask, newsize);
    negative_norm_scaled = imresize(neghist_norm, newsize).*mask_scaled;
    anypositive_norm_scaled = imresize(anyposhist_norm, newsize).*mask_scaled;
    negativeR_norm_scaled = imresize(neghistR_norm, [1 newsize(2)]);
    negativeZ_norm_scaled = imresize(neghistZ_norm, [newsize(1) 1]);
    r_scaled = imresize(Rgl(1,:), [1 newsize(2)]);
    z_scaled = imresize(Zgl(:,1), [newsize(1) 1]);
    for ci = 1:3
        positive_norm_scaled{ci} = imresize(poshist_norm{ci}, newsize).*mask_scaled;
        positiveR_norm_scaled{ci} = imresize(poshistR_norm{ci}, [1 newsize(2)]);
        positiveZ_norm_scaled{ci} = imresize(poshistZ_norm{ci}, [newsize(1) 1]);
    end
    
    % outline
    [B,~] = bwboundaries(mask_scaled,'noholes');
    
    % collect results
    result = struct('density', totaldensity,...
                    'density_scaled', density_scaled,...
                    'anypositive_norm_scaled', anypositive_norm_scaled,...
                    'negative_norm_scaled', negative_norm_scaled,...
                    'negativeR_norm_scaled', negativeR_norm_scaled,...
                    'negativeZ_norm_scaled', negativeZ_norm_scaled,...
                    'positive_norm_scaled', {positive_norm_scaled},...
                    'positiveR_norm_scaled', {positiveR_norm_scaled},...
                    'positiveZ_norm_scaled', {positiveZ_norm_scaled},...
                    'densitymask', mask,...
                    'densitymask_scaled',mask_scaled,...
                    'R',R,'Z',Z,...
                    'r_scaled',r_scaled,'z_scaled',z_scaled,...
                    'Rmax',Rmax,...
                    'Zmax',Zmax,...
                    'boundaries',{B});
end
