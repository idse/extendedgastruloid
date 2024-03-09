function res = plotRadialProfiles(stats, meta, options)
% average profile with equal number of datapoints per 'bin'

if ~isfield(options,'normalize')
    options.normalize = false;
end
if ~isfield(options,'nucChannels') && ~isfield(options,'cytChannels') && ~isfield(options,'ratioChannels')
    error('please specify channels');
end
if ~isfield(options,'nucChannels')
    options.nucChannels = [];
end
if ~isfield(options,'cytChannels')
    options.cytChannels = [];
end
if ~isfield(options,'ratioChannels')
    options.ratioChannels = [];
end
if ~isfield(options,'legend')
    options.legend = true;
end
if ~isfield(options,'colors')
    options.colors = lines(5);
    options.colors = options.colors([2 5 1 3 4],:);
end
if ~isfield(options,'std') %|| options.normalize == true
    options.std = false;
end
if ~isfield(options,'stdMethod')
    options.stdMethod = 'perColony';
end
if ~isfield(options,'Ilimmargin')
    options.Ilimmargin = 3;
end
if isfield(options,'FontSize')
    fs = options.FontSize;
else
    fs = 20;
end
if isfield(options,'conditionIdx')
    condi = options.conditionIdx;
else
    error('please specify options.conditionIdx');
end
if isfield(options,'pointsPerBin')
    ptsperbin = options.pointsPerBin;
else
    ptsperbin = 400;
end
if ~isfield(options,'useTrueRadius')
    options.useTrueRadius = true;
end
if ~isfield(options,'edgeDistance')
    options.edgeDistance = true;
end
if ~isfield(options,'time')
    time = size(stats.XY,2);
else
    time = options.time;
end
if ~isfield(options,'ratioMode')
    options.ratioMode = 'N:C';
end
if ~isfield(options,'interpmethod')
    interpmethod = 'linear';
else
    interpmethod = options.interpmethod;
end
if isfield(options, 'log1p')
    log1p = options.log1p;
else
    log1p = false;
end

nc = meta.nChannels;

% average with equal number of datapoints per 'bin' for all colonies
% combined
% this reduces artefacts due to small numbers of junk points on the
% colony edge
R = sqrt(sum(stats.XY{condi,time}.^2,2))*meta.xres;
nL = stats.nucLevel{condi,time};
cL = stats.cytLevel{condi,time};
if strcmp(options.ratioMode,'N:C')
    ratio = nL./cL; %nuc-cyt ratio; nL and cL are already background subtracted
    rpref = "N:C ratio ";
    disp(rpref);
elseif strcmp(options.ratioMode,'C:N')
    ratio = cL./nL; %nuc-cyt ratio; nL and cL are already background subtracted
    rpref = "C:N ratio ";
    disp(rpref);
end

% we usually normalize to threshold 
% and if log1p, then log(1+threshold) = 1
% so threshold = exp(1)-1 
if log1p
    for ci = 1:4
        %norm = log2(stats.thresholds(ci))/(exp(1)-1);     
        %nL(:,ci) = log(1+nL(:,ci)./norm);
        nL(:,ci) = log2(nL(:,ci));
    end
end

% set up a regularly spaced grid for interpolation
dr = 3; % step in micron
nominalradius = stats.radiusMicron{condi};
margin = 50;
maxradius = nominalradius + margin;
Ngrid = round(maxradius/dr);
radialgrid = linspace(0, maxradius, Ngrid);

% calculate profiles of individual colonies
%------------------------------------------------------------------
nuc_profiles = {};
nuc_profile_stds = {};
cyt_profiles = {};
cyt_profile_stds = {};
ratio_profiles = {};
ratio_profile_stds = {};
radii = {};

samples = unique(stats.sample{condi})';
for si = samples

    sbi = stats.sample{condi,time} == si;
    Rs = R(sbi);
    nLs = nL(sbi,:);
    cLs = cL(sbi,:);
    ratios = ratio(sbi,:);

    % estimate cell density and colony radius (which may differ
    % slightly from nominal radius of micropattern)
    PI = 3.1415;

    bw = 5;
    [f,xi] = ksdensity(Rs,'Bandwidth',bw,'BoundaryCorrection','reflection');

    dx = xi(2)-xi(1);
    Ncells = numel(Rs);
    celldensity{si} = Ncells*f./(2*PI*xi);

    nominalradius = stats.radiusMicron{condi};
    xrangeinside = xi > bw & xi < 0.9*nominalradius;
    avgdensityinside = sum(2*PI*celldensity{si}(xrangeinside).*xi(xrangeinside)*dx)/sum(2*PI*xi(xrangeinside)*dx);
    % define the radius as the point where the density drops to 10% of
    % the mean inside
    fullxrange = xi > bw & celldensity{si} > 0.1*avgdensityinside;
    trueradius = max(xi(fullxrange));

    if options.useTrueRadius
        radius = trueradius;
        if time == 1
            disp('adjusting for true radius of colony');
        end
    else
        radius = nominalradius;
    end
    radii{si} = radius;

    % calculate radial profiles of markers for individual colonies
    ptsperbinindividualcolonies = 10;
    [~,I] = sort(Rs);
    N = round(numel(I)/ptsperbinindividualcolonies);
    edges = round(linspace(1, numel(I), N));
    overlap = ptsperbinindividualcolonies; % overlap bins for smoothness

    r_tmp = zeros([N-1 1]);
    nuc_profile_tmp = zeros([N-1 nc]);
    nuc_profile_std_tmp = zeros([N-1 nc]);
    cyt_profile_tmp = zeros([N-1 nc]);
    cyt_profile_std_tmp = zeros([N-1 nc]);
    ratio_profile_tmp = zeros([N-1 nc]);
    ratio_profile_std_tmp = zeros([N-1 nc]);

    for i = 1:N-1
        ptidx = max(min(edges),edges(i)-overlap):min(edges(i+1)+overlap,max(edges));
        r_tmp(i) = mean(Rs(I(ptidx)));
        nuc_profile_tmp(i,:) = nanmean(nLs(I(ptidx),:));
        nuc_profile_std_tmp(i,:) = nanstd(nLs(I(ptidx),:));%/sqrt(ptsperbin);
        cyt_profile_tmp(i,:) = nanmean(cLs(I(ptidx),:));
        cyt_profile_std_tmp(i,:) = nanstd(cLs(I(ptidx),:));%/sqrt(ptsperbin);
        ratio_profile_tmp(i,:) = nanmean(ratios(I(ptidx),:));
        ratio_profile_std_tmp(i,:) = nanstd(ratios(I(ptidx),:));%/sqrt(ptsperbin);
    end

    % convert radius to edge distance
    if options.edgeDistance == true
        r_tmp = radius - r_tmp;
        R(sbi) = radius - R(sbi);
    end

    nuc_profiles{si} = zeros([Ngrid nc]);
    nuc_profile_stds{si} = zeros([Ngrid nc]);
    cyt_profiles{si} = zeros([Ngrid nc]);
    cyt_profile_stds{si} = zeros([Ngrid nc]);
    ratio_profiles{si} = zeros([Ngrid nc]);
    ratio_profile_stds{si} = zeros([Ngrid nc]);

    % interpolate everything on a fixed radial grid
    if N > 2 %skip empty or nearly empty colonies
        for ci = 1:nc
            nuc_profiles{si}(:,ci) = interp1(r_tmp, nuc_profile_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            cyt_profiles{si}(:,ci) = interp1(r_tmp, cyt_profile_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            nuc_profile_stds{si}(:,ci) = interp1(r_tmp, nuc_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            cyt_profile_stds{si}(:,ci) = interp1(r_tmp, cyt_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            ratio_profiles{si}(:,ci) = interp1(r_tmp, ratio_profile_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
            ratio_profile_stds{si}(:,ci) = interp1(r_tmp, ratio_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
        end
        celldensity{si} = interp1(radius-xi, celldensity{si}, radialgrid, interpmethod,'extrap');
    else
        nuc_profiles{si} = NaN([Ngrid nc]);
        cyt_profiles{si} = NaN([Ngrid nc]);
        nuc_profile_stds{si} = NaN([Ngrid nc]);
        cyt_profile_stds{si} = NaN([Ngrid nc]);
        ratio_profiles{si} = NaN([Ngrid nc]);
        ratio_profile_stds{si} = NaN([Ngrid nc]);
        celldensity{si} = NaN(1,Ngrid);
    end
end

% average over all colonies with equal number of datapoints per 'bin'
%------------------------------------------------------------------
[~,I] = sort(R);

N = round(numel(I)/ptsperbin);
edges = round(linspace(1, numel(I), N));
overlap = round(ptsperbin/2);

r_tmp = zeros([N-1 1]);
nuc_profile_tmp = zeros([N-1 nc]);
nuc_profile_std_tmp = zeros([N-1 nc]);
cyt_profile_tmp = zeros([N-1 nc]);
cyt_profile_std_tmp = zeros([N-1 nc]);
%also do N:C ratio
ratio_profile_tmp = zeros(N-1,nc);
ratio_profile_std_tmp = zeros(N-1,nc);

for i = 1:N-1
    ptidx = max(min(edges),edges(i)-overlap):min(edges(i+1)+overlap,max(edges));
    r_tmp(i) = mean(R(I(ptidx)));
    nuc_profile_tmp(i,:) = nanmean(nL(I(ptidx),:));
    nuc_profile_std_tmp(i,:) = nanstd(nL(I(ptidx),:));%/sqrt(ptsperbin);
    cyt_profile_tmp(i,:) = nanmean(cL(I(ptidx),:));
    cyt_profile_std_tmp(i,:) = nanstd(cL(I(ptidx),:));%/sqrt(ptsperbin);
    %N:C ratio
    ratio_profile_tmp(i,:) = nanmean(ratio(I(ptidx),:));
    ratio_profile_std_tmp(i,:) = nanstd(ratio(I(ptidx),:));
end

% interpolate on grid
nuc_profile = zeros([Ngrid nc]);
nuc_profile_std = zeros([Ngrid nc]);
cyt_profile = zeros([Ngrid nc]);
cyt_profile_std = zeros([Ngrid nc]);
%N:C ratio
ratio_profile = zeros([Ngrid nc]);
ratio_profile_std = zeros([Ngrid nc]);
for ci = 1:nc
    nuc_profile(:,ci) = interp1(r_tmp, nuc_profile_tmp(:,ci)', radialgrid, interpmethod, 'extrap')';
    cyt_profile(:,ci) = interp1(r_tmp, cyt_profile_tmp(:,ci)', radialgrid, interpmethod,'extrap')';
    nuc_profile_std(:,ci) = interp1(r_tmp, nuc_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
    cyt_profile_std(:,ci) = interp1(r_tmp, cyt_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
    %N:C ratio
    ratio_profile(:,ci) = interp1(r_tmp, ratio_profile_tmp(:,ci)', radialgrid, interpmethod,'extrap')';
    ratio_profile_std(:,ci) = interp1(r_tmp, ratio_profile_std_tmp(:,ci)', radialgrid,interpmethod,'extrap')';
end

% alternative error bar: std of colonies
% since std from 4 colonies is very noisy but true std should vary
% smoothly in neighboring points we average the estimate of the std
% between neighbors
nuc_profile_colstd = imfilter(nanstd(cat(3, nuc_profiles{:}),[],3),[1 1 1 1 1]'/5, 'replicate');
cyt_profile_colstd = imfilter(nanstd(cat(3, cyt_profiles{:}),[],3),[1 1 1 1 1]'/5, 'replicate');
ratio_profile_colstd = imfilter(nanstd(cat(3, ratio_profiles{:}),[],3),[1 1 1 1 1]'/5, 'replicate');

% cut off the plot where there is no data
maxtrueradius = max([radii{:}]);
if options.edgeDistance == false
    nuc_profile(radialgrid > maxtrueradius, :) = NaN;
    cyt_profile(radialgrid > maxtrueradius, :) = NaN;
    ratio_profile(radialgrid > maxtrueradius, :) = NaN;
    nuc_profile_colstd(radialgrid > maxtrueradius, :) = NaN;
    cyt_profile_colstd(radialgrid > maxtrueradius, :) = NaN;
    ratio_profile_colstd(radialgrid > maxtrueradius, :) = NaN;
end

% determine limits for normalization
if isfield(options,'nuclimits')
    nuclimits = options.nuclimits;
else
    nuclimits = [min(nuc_profile(1:end,:))' max(nuc_profile(1:end,:))'];
end
if isfield(options,'cytlimits')
    cytlimits = options.cytlimits;
else
    cytlimits = [min(cyt_profile)' max(cyt_profile)'];
end
if isfield(options,'ratiolimits')
    ratiolimits = options.ratiolimits;
else
    ratiolimits = [min(ratio_profile(1:end,:))' max(ratio_profile(1:end,:))'];
end

if options.normalize
    for ci = 1:size(nuc_profile,2)
        nuc_profile(:,ci) = (nuc_profile(:,ci) - nuclimits(ci,1))./(nuclimits(ci,2) - nuclimits(ci,1));
        cyt_profile(:,ci) = (cyt_profile(:,ci) - cytlimits(ci,1))./(cytlimits(ci,2) - cytlimits(ci,1));
        nuc_profile_std(:,ci) = nuc_profile_std(:,ci)./(nuclimits(ci,2) - nuclimits(ci,1));
        cyt_profile_std(:,ci) = cyt_profile_std(:,ci)./(cytlimits(ci,2) - cytlimits(ci,1));
        nuc_profile_colstd(:,ci) = nuc_profile_colstd(:,ci)./(nuclimits(ci,2) - nuclimits(ci,1));
        cyt_profile_colstd(:,ci) = cyt_profile_colstd(:,ci)./(cytlimits(ci,2) - cytlimits(ci,1));
        %N:C ratio
        ratio_profile(:,ci) = (ratio_profile(:,ci) - ratiolimits(ci,1))./(ratiolimits(ci,2) - ratiolimits(ci,1));
        ratio_profile_std(:,ci) = ratio_profile_std(:,ci)./(ratiolimits(ci,2) - ratiolimits(ci,1));
        ratio_profile_colstd(:,ci) = ratio_profile_colstd(:,ci)./(ratiolimits(ci,2) - ratiolimits(ci,1));
    end
end

r = radialgrid;
legendentries = {};
hold on
if ~isempty(options.nucChannels)
    for cii = 1:numel(options.nucChannels)
        ci = options.nucChannels(cii)+1;
        plot(r, nuc_profile(:,ci),'LineWidth',3, 'Color', options.colors(cii,:))

        legendentries = [legendentries, ['nuc ' meta.channelLabel{ci}],[]];
    end
    for cii = 1:numel(options.nucChannels)
        ci = options.nucChannels(cii)+1;
        if options.std
            if strcmp(options.stdMethod,'perColony')
                %errorbar(r, nuc_profile(:,ci), nuc_profile_colstd(:,ci),'LineWidth',1, 'Color', options.colors(cii,:),'HandleVisibility','off');
                err = nuc_profile_colstd(:,ci);
            elseif strcmp(options.stdMethod,'neighborCells')
                err = nuc_profile_std(:,ci);
                %errorbar(r, nuc_profile(:,ci), nuc_profile_std(:,ci),'LineWidth',1, 'Color', options.colors(cii,:),'HandleVisibility','off');
            end
            
            plerr = nuc_profile(:,ci) + err;
            minerr = nuc_profile(:,ci) - err;
            good = ~isnan(plerr);
            fill([r(good),fliplr(r(good))],[plerr(good)', fliplr(minerr(good)')],...
                options.colors(cii,:),'FaceAlpha',0.2,'EdgeColor','none');
        end
    end
    if options.log1p
        ylabelstr = 'log( intensity )';
    else
        ylabelstr = 'intensity';
    end
    if options.normalize
        ylabelstr = [ylabelstr ' (a.u.)'];
    end
    ylabel(ylabelstr);
end

if ~isempty(options.cytChannels)
    for cii = numel(options.cytChannels):-1:1
        ci = options.cytChannels(cii) + 1;

        plot(r, cyt_profile(:,ci),'LineWidth',3, 'Color', options.colors(ci,:))

        legendentries = [legendentries, ['cyt ' meta.channelLabel{ci}]];
    end
    for cii = numel(options.cytChannels):-1:1
        ci = options.cytChannels(cii) + 1;
        if options.std
            if strcmp(options.stdMethod,'perColony')
                errorbar(r, cyt_profile(:,ci), cyt_profile_colstd(:,ci),'LineWidth',1, 'Color', options.colors(ci,:),'HandleVisibility','off');
            elseif strcmp(options.stdMethod,'neighborCells')
                errorbar(r, cyt_profile(:,ci), cyt_profile_std(:,ci),'LineWidth',1, 'Color', options.colors(ci,:),'HandleVisibility','off');
            end
        end
    end
    if options.log1p
        ylabelstr = 'log(intensity)';
    else
        ylabelstr = 'intensity';
    end
    if options.normalize
        ylabelstr = [ylabelstr ' (a.u.)'];
    end
    ylabel(ylabelstr);
end

if ~isempty(options.ratioChannels)
    for cii = numel(options.ratioChannels):-1:1
        ci = options.ratioChannels(cii) + 1;

        plot(r, ratio_profile(:,ci),'LineWidth',3, 'Color', options.colors(ci,:))

        legendentries = [legendentries, strcat(rpref, meta.channelLabel{ci})];
    end
    for cii = numel(options.ratioChannels):-1:1
        ci = options.ratioChannels(cii) + 1;
        if options.std
            if strcmp(options.stdMethod,'perColony')
                errorbar(r, ratio_profile(:,ci), ratio_profile_colstd(:,ci),'LineWidth',1, 'Color', options.colors(ci,:),'HandleVisibility','off');
            elseif strcmp(options.stdMethod,'neighborCells')
                errorbar(r, ratio_profile(:,ci), ratio_profile_std(:,ci),'LineWidth',1, 'Color', options.colors(ci,:),'HandleVisibility','off');
            end
        end
    end
    if options.normalize
        ylabel('ratio (a.u.)')
    else
        ylabel('ratio')
    end
end
hold off

if options.legend
    legend(legendentries, 'Location','best','AutoUpdate','off');
end

% make it pretty
fgc = 'k';
bgc = 'w';
graphbgc = 1*[1 1 1];

if options.edgeDistance
    xlabel('edge distance (um)', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
else
    xlabel('r (um)', 'FontSize',fs,'FontWeight','Bold','Color',fgc)
end

set(gcf,'color',bgc);
set(gca, 'LineWidth', 2);
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca,'XColor',fgc);
set(gca,'YColor',fgc);
set(gca,'Color',graphbgc);

datatable = table(r', nuc_profile, nuc_profile_std, nuc_profile_colstd, cyt_profile_colstd,...
            cyt_profile, cyt_profile_std,ratio_profile, ratio_profile_std,...
            'VariableNames', {'r', 'nuc_profile', 'nuc_profile_std',...
            'nuc_profile_colstd', 'cyt_profile_colstd',...
            'cyt_profile', 'cyt_profile_std',...
            'ratio_profile', 'ratio_profile_std'});

res = struct('r',r,'nuc_profile',nuc_profile, 'nuc_profile_std',nuc_profile_std,...
    'cyt_profile', cyt_profile, 'cyt_profile_std', cyt_profile_std,...
    'ratio_profile',ratio_profile, 'ratio_profile_std',ratio_profile_std,...
    'celldensity',{celldensity},'nuc_profile_colstd',nuc_profile_colstd,...
    'cyt_profile_colstd',cyt_profile_colstd,'trueradius',{radii}, 'datatable', datatable);
end

