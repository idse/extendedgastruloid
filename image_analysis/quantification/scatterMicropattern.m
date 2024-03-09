function results = scatterMicropattern(stats, meta, dataDir, options)

% options.channelMax in units of the norm

if ~isfield(options,'channelCombos')
    channelCombos = {[3 2 4], [4 3 2], [4 2 3]};
else
    channelCombos = options.channelCombos;
    if min(cellfun(@numel, options.channelCombos)) < 3
        error('channel combos needs 3 channels, e.g. {[3 2 4]}');
    end
end
channelThresholds = options.channelThresholds;

% options.norm: 
% norm : 1 x Nchannels x 2 or Nconditions x Nchannels x 2, if not specified norm=threshold
% last dimension 2 values for min and max vals
% if Nconditions x Nchannels x 2, use normIdx to specify each conditions to use
% its own norm or it will use the first

% potential different norm for different conditions (necessary when
% combining different datasets)
% vector for each condition in stats with condition from which to get the
% norm (based on thresholds)

if ~isfield(options,'normIdx')
    normIdx = zeros([meta.nWells 1],'uint16') + 1;
else
    normIdx = options.normIdx;
end

% log1p : boolean vector to transform each channel or not
if isfield(options, 'log1p')
    log1p = options.log1p;
else
    log1p = true([1 meta.nChannels]);
end

if isfield(options, 'edgeDistance')
    edgeDistance = options.edgeDistance;
else
    edgeDistance = true;
end

if isfield(options, 'conditionsCombined')
    conditionsCombined = options.conditionsCombined;
else
    conditionsCombined = true;
end

if isfield(options, 'selectDisType')
    DisType = options.selectDisType;
else
    DisType = [1, 1, 1, 1];
end
for i = 1:length(options.axisLabel)
    if DisType(i) == 2
        options.axisLabel{i} = ['cyto' options.axisLabel{i}];
    elseif DisType(i) == 3
        options.axisLabel{i} = ['NCratio' options.axisLabel{i}];
    elseif DisType(i) == 4
        options.axisLabel{i} = ['CNratio' options.axisLabel{i}];
    end
end
if ~isfield(options, 'showThreshold')
    options.showThreshold = true([1 meta.nChannels]);
end
if ~isfield(options, 'minimal')
    options.minimal = false;
end
if isfield(options, 'checkPos') 
    positive = options.checkPos;
else
    positive = cell(numel(meta.conditions),1);
    for i = 1:numel(meta.conditions)
        positive{i} = zeros(size(stats.nucLevel{i},1),1);
    end
end

% NOTE : FIRST BLOCK OF CODE HAS BEEN UPDATED, BUT CONDITIONS COMBINED
% PLOTS BELOW STILL USE OLDER CODE WHICH IS LESS FLEXIBLE (CAN'T TURN OF
% THRESHOLDS FOR SOME CHANNELS ETC)
% WE CAN FIX THAT WHEN WE NEED TO

% case 1: color for edge distance, 
% case 2: color for 3rd gene
% case 3: density of points (like FACS)

if isfield(options,'colorDensity') && options.colorDensity == false
    cases = [1 2];
else
    cases = [1 2 3];
end

if strcmp(options.imageType,'MP')
    cases = [1 2 3];
elseif strcmp(options.imageType,'Disordered')
    cases = [1 3];
end

for casei = cases
for conditionIdx = options.conditionIdx

    % see definition of normidx above
    ni = normIdx(conditionIdx);

    if casei == 2
        XY = stats.XY{conditionIdx}(:,1:2);
        if edgeDistance
            dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;
            dist(dist<0) = 0; % cut off what is outside the true radius (shouldn't be much)
        else
            dist = sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;
            dist(dist > options.radiusMicron) = options.radiusMicron;
        end
    end
    
    for comboi = 1:numel(channelCombos)
        
        figure('Position',[0 0 700 500])
        
        channelIdx = channelCombos{comboi};

        xmin = [0 0 0];
        xmax = [0 0 0]; % xmax = [xmax ymax zmax], etc
        xthresh = [0 0 0];
        norm = [0 0; 0 0; 0 0]; % first column: min, second: max
        ticklist = {};
        x = {};
        % for each channel we are graphing
        for ci = 1:3
            
            % set norm
            % we usually normalize to threshold 
            % and if log1p, then log(1+threshold) = 1
            % so threshold = exp(1)-1 
            if log1p(channelIdx(ci))
                if isfield(options,'norm')
                    norm(ci,1) = options.norm(ni,channelIdx(ci),1)/(exp(1)-1);
                    norm(ci,2) = options.norm(ni,channelIdx(ci),2)/(exp(1)-1);
                else
                    norm(ci,2) = options.channelThresholds(ni,channelIdx(ci))/(exp(1)-1);     
                end
                ticklist{ci} = [0 1 2 3 4];
            else
                if isfield(options,'norm')
                    norm(ci,1) = options.norm(ni, channelIdx(ci),1);
                    norm(ci,2) = options.norm(ni, channelIdx(ci),2);
                else
                    norm(ci,2) = options.channelThresholds(ni,channelIdx(ci));     
                end
                ticklist{ci} = [0 5 10 20 30 40 50];
            end
            
            % define each of the coordinates
            if DisType(channelIdx(ci)) == 1
                x{ci} = (stats.nucLevel{conditionIdx}(:,channelIdx(ci))-norm(ci,1))/norm(ci,2);
                xP{ci} = (stats.nucLevel{conditionIdx}(:,channelIdx(ci))-norm(ci,1)).*positive{conditionIdx}/norm(ci,2);
            elseif DisType(channelIdx(ci)) == 2
                x{ci} = (stats.cytLevel{conditionIdx}(:,channelIdx(ci))-norm(ci,1))/norm(ci,2);
                xP{ci} = (stats.cytLevel{conditionIdx}(:,channelIdx(ci))-norm(ci,1)).*positive{conditionIdx}/norm(ci,2);
            elseif DisType(channelIdx(ci)) == 3
                x{ci} = ((stats.nucLevel{conditionIdx}(:,channelIdx(ci))./stats.cytLevel{conditionIdx}(:,channelIdx(ci)))-norm(ci,1))/norm(ci,2);
                xP{ci} = ((stats.nucLevel{conditionIdx}(:,channelIdx(ci))./stats.cytLevel{conditionIdx}(:,channelIdx(ci)))-norm(ci,1)).*positive{conditionIdx}/norm(ci,2);
            elseif DisType(channelIdx(ci)) == 4
                x{ci} = ((stats.cytLevel{conditionIdx}(:,channelIdx(ci))./stats.nucLevel{conditionIdx}(:,channelIdx(ci)))-norm(ci,1))/norm(ci,2);
                xP{ci} = ((stats.cytLevel{conditionIdx}(:,channelIdx(ci))./stats.nucLevel{conditionIdx}(:,channelIdx(ci)))-norm(ci,1)).*positive{conditionIdx}/norm(ci,2);
            end
            % limits on the graph
            %xmin(ci) = (stats.lim{channelIdx(ci)}(1)-norm(ci,1))/norm(ci,2);
            if isfield(options, 'channelMax')
                xmax(ci) = options.channelMax(channelIdx(ci));
            else
                xmax(ci) = (stats.lim{channelIdx(ci)}(2)-norm(ci,1))/norm(ci,2);
            end
            
            % thresholds
            % if there are different threshold for each condition (should
            % typically not be used)
            if size(channelThresholds,1) > 1
                xthresh(ci) = (channelThresholds(conditionIdx,channelIdx(ci))-norm(ci,1))/norm(ci,2);
            else
                xthresh(ci) = (channelThresholds(channelIdx(ci))-norm(ci,1))/norm(ci,2);
            end
            
            % log transform
            if log1p(channelIdx(ci))
                x{ci} = log(1+x{ci});
                xP{ci} = log(1+xP{ci});
                xmin(ci) = log(1+xmin(ci));
                xmax(ci) = log(1+xmax(ci));
            	xthresh(ci) = log(1+xthresh(ci));
            end
        end
        
        % for 3rd coordinate set cutoff (on colorbar)
        x{3}(x{3}>xmax(3)) = xmax(3);
        if ~log1p
            x{3}(x{3}<0) = 0;
        end
        
        if casei == 1
            c = x{3};
            %c = stats.nucLevel{conditionIdx}(:,cidx3) > options.channelThresholds(ni,cidx3);
            suffix = ['scatter_' meta.channelLabel{channelIdx(3)}];
            if ~isfield(options,'axisLabel')
                collabel = meta.channelLabel{channelIdx(3)};
            else
                collabel = options.axisLabel{channelIdx(3)};
            end
            
        elseif casei == 2
            c = dist;
            suffix = 'scatter_dist';
            if edgeDistance
                collabel = 'distance from edge';
            else
                collabel = 'r (\mu m)';
            end
            
        elseif casei == 3
            x{1} = real(x{1});
            x{2} = real(x{2});
            [f,xi] = ksdensity([x{1},x{2}]);
            np = sqrt(size(xi,1));
            Xi = reshape(xi(:,1),[np,np]); Yi = reshape(xi(:,2),[np,np]); F = reshape(f,[np,np]);
            z = interp2(Xi,Yi,F,x{1},x{2});

            tol = 0.02;
            n = length(z);
            zs = sort(z);
            zmin = zs(max(ceil(n*tol),1));
            zmax = zs(floor(n*(1-tol)));
            z(z < zmin) = zmin; z(z > zmax) = zmax;
            density = z;
            c = z;
            collabel = 'cell density';
            suffix = 'scatter_density';
        end
        
        if ~options.minimal
            lw = 3;
            fs = 20;
            pfs = 20;
            s = 10;
        else
            lw = 5;
            fs = 50;
            pfs = 20;
            s = 20;
        end
        
        h = scatter(x{1},x{2}, s, c,'filled');
        hold on
        if isfield(options, 'checkPos')
        hh = scatter(xP{1},xP{2}, s,'filled', 'magenta');
        hold on
        end
        
        try
            colormap(turbo)
            %colormap(brewermap([],"Reds"))
            %colormap(flipud(brewermap([],"RdBu")))
            %colormap(flipud(brewermap([],"RdGy")))
        catch
            colormap(jet)
        end
        
        if ~options.minimal 
            c = colorbar;
            c.Label.String = collabel;
            if casei == 3 % the numbers for density have no meaning
                c.TickLabels = [];
            end
        end
        
        xlim([xmin(1) xmax(1)])
        ylim([xmin(2) xmax(2)])
        
        %xlim([min(x{1}), max(x{1})])
        %ylim([min(x{2}), max(x{2})])

        if ~isfield(options,'axisLabel')
            xlabel([meta.channelLabel{channelIdx(1)}])
            ylabel([meta.channelLabel{channelIdx(2)}])
        else
            xlabel([options.axisLabel{channelIdx(1)}])
            ylabel([options.axisLabel{channelIdx(2)}])
        end
        set(gca, 'LineWidth', lw);
        set(gca,'FontSize', fs)
        set(gca,'FontWeight', 'bold')
        set(gca,'Color','w');
        axis square;

        
        if options.minimal
            xticks([]);
            yticks([]);
        else
            xticks(ticklist{1});
            yticks(ticklist{2});
        end
        %ax = ancestor(h, 'axes');
        %ax.XAxis.Exponent = 3;
        %ax.YAxis.Exponent = 3;

        % +-, ++, -+, --
        pnidx = x{2} > xthresh(2) & x{1} <= xthresh(1);
        ppidx = x{2} > xthresh(2) & x{1} > xthresh(1);
        npidx = x{2} <= xthresh(2) & x{1} > xthresh(1);
        nnidx = x{2} <= xthresh(2) & x{1} <= xthresh(1);
        pn = 100*sum(pnidx)/numel(x{2});
        pp = 100*sum(ppidx)/numel(x{2});
        np = 100*sum(npidx)/numel(x{2});
        nn = 100*sum(nnidx)/numel(x{2});
        
        disp(['y+x-: ' num2str(pn,3) '%']);
        disp(['y+x+: ' num2str(pp,3) '%']);
        disp(['y-x+: ' num2str(np,3) '%']);
        disp(['y-x-: ' num2str(nn,3) '%']);
        fc = 'k';

        if options.showThreshold(channelIdx(1)) 
            line(xthresh(1)*[1 1],[xmin(2) xmax(2)],'LineStyle','--','LineWidth',lw,'Color','k');
        end
        if options.showThreshold(channelIdx(2)) 
            line([xmin(1) xmax(1)],xthresh(2)*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
        end
        if ~options.minimal
            title(meta.conditions{conditionIdx});
        end

        if ~options.minimal 
            
            if options.showThreshold(channelIdx(1)) && options.showThreshold(channelIdx(2)) 
                text(xmin(1) + (xthresh(1) - xmin(1))/2, xmin(2) + (xthresh(2)-xmin(2))/2, [num2str(nn,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
                text(xmin(1) + (xthresh(1) - xmin(1))/2, xthresh(2) + (xmax(2) - xthresh(2))/2, [num2str(pn,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
                text(xthresh(1) + (xmax(1) - xthresh(1))/2, xmin(2) + (xthresh(2)-xmin(2))/2, [num2str(np,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
                text(xthresh(1) + (xmax(1) - xthresh(1))/2, xthresh(2) + (xmax(2) - xthresh(2))/2, [num2str(pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
            elseif options.showThreshold(channelIdx(1)) && ~options.showThreshold(channelIdx(2)) 
                text(xmin(1) + (xthresh(1)-xmin(1))/2, xthresh(2), [num2str(nn + pn,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
                text(xthresh(1) + (xmax(1) - xthresh)/2, xthresh(2), [num2str(np + pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
            elseif options.showThreshold(channelIdx(2)) && ~options.showThreshold(channelIdx(1)) 
                text(xthresh(1), xthresh(2)/2, [num2str(nn + np,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
                text(xthresh(1), xthresh(2)  + (xmax(2) - xthresh(2))/2, [num2str(pn + pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
            end
            
            
            N = size(stats.nucLevel{conditionIdx},1);
            text(xmin(1) + 0.05*(xmax(1)-xmin(1)), xmax(2)*0.9, ['N=' num2str(N)],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
        end    
        saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' meta.conditions{conditionIdx} '_' suffix '.png'])); 
        close;
    end
end
end

% return results table to write to CSV
results = struct();
results.dataTable = table(x{1},x{2},x{3}, dist, density,'VariableNames', [meta.channelLabel(channelIdx), 'distance', 'density']);

%
if conditionsCombined

%--------------------------------------------------------------------------
% combine conditions in one plot and color for condition
%--------------------------------------------------------------------------

colors = lines(7);
colors = colors([1:3 5 4 6 7],:);
%colors = lines(numel(options.conditionIdx));

for comboi = 1:numel(channelCombos)
    
    figure('Position',[0 0 600 500])
    hold on;
    xall = [];
    yall = [];
    call = [];
    
    for i = 1:numel(options.conditionIdx)
        
        conditionIdx = options.conditionIdx(end+1-i);
        
        ni = normIdx(conditionIdx);
        posidx = stats.nucLevel{ni} < options.channelThresholds(ni,:);
        %norm = sum(stats.nucLevel{ni}.*posidx,1)./sum(posidx,1);
        norm = options.channelThresholds(ni,:)/(exp(1)-1);
        
        channelIdx = channelCombos{comboi}(1:2);
        
        if isfield(options, 'channelMax')
            xmax = options.channelMax(channelIdx(1));
            ymax = options.channelMax(channelIdx(2));
        else
            xmax = stats.lim{channelIdx(1)}(2)/norm(channelIdx(1));
            ymax = stats.lim{channelIdx(2)}(2)/norm(channelIdx(2));
        end
        
        if size(channelThresholds,1) > 1
            xthresh = channelThresholds(conditionIdx,channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(conditionIdx,channelIdx(2))/norm(channelIdx(2));
        else
            xthresh = channelThresholds(channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(channelIdx(2))/norm(channelIdx(2));
        end

        XY = stats.XY{conditionIdx};
        dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;

        %hist(dist,50)

        x = stats.nucLevel{conditionIdx}(:,channelIdx(1))/norm(channelIdx(1));
        y = stats.nucLevel{conditionIdx}(:,channelIdx(2))/norm(channelIdx(2));
        
        if numel(channelCombos{comboi}) == 3
            cidx3 = channelCombos{comboi}(3);
        else
            cidx3 = setdiff(2:4,channelIdx);
            % below is in case one of the indices is DAPI so that there are 
            % multiple channels remaining, which is usually not the case
            cidx3 = cidx3(1); 
        end
        
        z = stats.nucLevel{conditionIdx}(:,cidx3)/norm(cidx3);
        %z = stats.nucLevel{conditionIdx}(:,4);
        c = repmat(colors(i,:),[size(x,1) 1]);
        %cmax = 200;
        %c(c > cmax) = cmax;
        %c = stats.sample{conditionIdx};
        
        % this is for the code below to randomize points from different
        % sets (not used right now)
        xall = cat(1,xall,x);
        yall = cat(1,yall,y);
        call = cat(1,call,c);
        
        xlims = stats.lim{channelIdx(1)}/norm(channelIdx(1));
        ylims = stats.lim{channelIdx(2)}/norm(channelIdx(2));
        xlims(2) = xmax;
        ylims(2) = ymax;
        
        if log1p
            x = log(1+x); 
            y = log(1+y);
            xmax = log(1+xmax);
            ymax = log(1+ymax);
            xthresh = log(1+xthresh);
            ythresh = log(1+ythresh);
            xlims = log(1+xlims);
            ylims = log(1+ylims);
        end

        ma = 1;
        h = scatter(x,y, s, c,'filled','MarkerFaceAlpha',ma);
    end
    
%     shuffleidx = randperm(size(xall,1));
%     xall = xall(shuffleidx);
%     yall = yall(shuffleidx);
%     call = call(shuffleidx,:);
%     
%     if log1p
%         xall = log(1+xall); 
%         yall = log(1+yall);
%     end
%     h = scatter(xall,yall, s, call,'filled','MarkerFaceAlpha',0.5);

    colormap(turbo)
    %c = colorbar;
    %c.Label.String = 'distance from edge';

    xlim(xlims)
    ylim(ylims)

    if ~isfield(options,'axisLabel')
        xlabel([meta.channelLabel{channelIdx(1)}])
        ylabel([meta.channelLabel{channelIdx(2)}])
    else
        xlabel([options.axisLabel{channelIdx(1)}])
        ylabel([options.axisLabel{channelIdx(2)}])
    end
    
    set(gca, 'LineWidth', lw);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'Color','w');
    axis square;
        
%     fc = 'k';
%         text(xthresh/2, ythresh/2, [num2str(mm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w');
%         text(xthresh/2, ythresh + (ymax - ythresh)/2, [num2str(pm,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
%         text(xthresh + (xmax - xthresh)/2, ythresh/2, [num2str(mp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
%         text(xthresh + (xmax - xthresh)/2, ythresh + (ymax - ythresh)/2, [num2str(pp,3) '%'],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')

%         N = size(stats.nucLevel{conditionIdx},1);
%         text(0, ymax*0.9, ['N=' num2str(N)],'FontSize',pfs,'Color',fc,'FontWeight', 'bold','Background','w')
    
    if options.showThreshold(channelIdx(1)) 
        line(xthresh*[1 1],[0 ymax],'LineStyle','--','LineWidth',lw,'Color','k');
    end
    if options.showThreshold(channelIdx(2)) 
        line([0 xmax],ythresh*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
    end
    %title(meta.conditions{conditionIdx});

    lfs = 20;
    if isfield(options,'legendstr')
        legend(options.legendstr,'FontSize',lfs,'Location','NorthWest');
    else
        legend(meta.conditions(options.conditionIdx(end:-1:1)),'FontSize',lfs,'Location','NorthWest');
    end
    hold off
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' num2str(options.conditionIdx) '_scatter_conditions.png'])); 
    close;
end

%------------------------------------------------------
% combine conditions in one plot and color for remaining marker or
% condition with shuffled points
%------------------------------------------------------

colors = lines(numel(options.conditionIdx));

for casei = cases

for comboi = 1:numel(channelCombos)
    
    figure('Position',[0 0 600 500])
    hold on;
    xall = [];
    yall = [];
    call = [];
    
    for i = 1:numel(options.conditionIdx)
        
        conditionIdx = options.conditionIdx(i);%end+1-i
        
        ni = normIdx(conditionIdx);
        %posidx = stats.nucLevel{ni} < options.channelThresholds(ni,:);
        %norm = sum(stats.nucLevel{ni}.*posidx,1)./sum(posidx,1);
        norm = options.channelThresholds(ni,:)/(exp(1)-1);
        
        channelIdx = channelCombos{comboi}(1:2);
        if isfield(options, 'channelMax')
            xmax = options.channelMax(channelIdx(1));
            ymax = options.channelMax(channelIdx(2));
        else
            xmax = stats.lim{channelIdx(1)}(2)/norm(channelIdx(1));
            ymax = stats.lim{channelIdx(2)}(2)/norm(channelIdx(2));
        end
        if size(channelThresholds,1) > 1
            xthresh = channelThresholds(conditionIdx,channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(conditionIdx,channelIdx(2))/norm(channelIdx(2));
        else
            xthresh = channelThresholds(channelIdx(1))/norm(channelIdx(1));
            ythresh = channelThresholds(channelIdx(2))/norm(channelIdx(2));
        end

        XY = stats.XY{conditionIdx};
        dist = options.radiusMicron - sqrt(pdist2(XY,[0 0],'squaredeuclidean'))*meta.xres;

        %hist(dist,50)

        x = stats.nucLevel{conditionIdx}(:,channelIdx(1))/norm(channelIdx(1));
        y = stats.nucLevel{conditionIdx}(:,channelIdx(2))/norm(channelIdx(2));
        
        if numel(channelCombos{comboi}) == 3
            cidx3 = channelCombos{comboi}(3);
        else
            cidx3 = setdiff(2:4,channelIdx);
            % below is in case one of the indices is DAPI so that there are 
            % multiple channels remaining, which is usually not the case
            cidx3 = cidx3(1); 
        end
        
        z = stats.nucLevel{conditionIdx}(:,cidx3)/norm(cidx3);
        zcut = 3;
        z(z<0) = 0;
        z(z>zcut) = zcut;
        
        %zscaled = mat2gray(z);
        %zcol = imadjust(zscaled, stretchlim(zscaled, [0 0.999]));
        
        if casei == 1
            c = z;
            suffix = 'scatter3_conditions';
        elseif casei == 2
            c = repmat(colors(i,:),[size(x,1) 1]);
            disp([meta.conditions{conditionIdx} num2str(colors(i,:))]);
            suffix = 'scatter_conditions_shuffle';
        end
        %cmax = 200;
        %c(c > cmax) = cmax;
        %c = stats.sample{conditionIdx};
        
        % this is for the code below to randomize points from different
        % sets 
        xall = cat(1,xall,x);
        yall = cat(1,yall,y);
        call = cat(1,call,c);
        
        xlims = stats.lim{channelIdx(1)}/norm(channelIdx(1));
        ylims = stats.lim{channelIdx(2)}/norm(channelIdx(2));
        xlims(2) = xmax;
        ylims(2) = ymax;
    end
    
    shuffleidx = randperm(size(xall,1));
    xall = xall(shuffleidx);
    yall = yall(shuffleidx);
    call = call(shuffleidx,:);
    
    if log1p
        xall = log(1+xall); 
        yall = log(1+yall);
        xmax = log(1+xmax);
        ymax = log(1+ymax);
        xthresh = log(1+xthresh);
        ythresh = log(1+ythresh);
        xlims = log(1+xlims);
        ylims = log(1+ylims);
    end
  
    h = scatter(xall,yall, s, call,'filled','MarkerFaceAlpha',ma);

    xlim(xlims)
    ylim(ylims)

    xlabel([meta.channelLabel{channelIdx(1)}])% ' (xm)'])
    ylabel([meta.channelLabel{channelIdx(2)}])% ' (xm)'])
    
    set(gca, 'LineWidth', lw);
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca,'Color','w');
    axis square;

    if casei == 1
        colormap(turbo)
        c = colorbar;
        c.Label.String = meta.channelLabel{cidx3};
    end
    %title(meta.conditions{conditionIdx});

    line(xthresh*[1 1],[0 ymax],'LineStyle','--','LineWidth',lw,'Color','k');
    line([0 xmax],ythresh*[1 1],'LineStyle','--','LineWidth',lw,'Color','k');
    
    hold off
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_' meta.channelLabel{channelIdx(2)} '_' num2str(options.conditionIdx) '_' suffix '.png'])); 
    close;
end
end
end

end