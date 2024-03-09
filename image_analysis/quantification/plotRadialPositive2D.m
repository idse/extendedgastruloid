function radpos2Dcomb = plotRadialPositive2D(positions, meta, channelThresholds, options)

    if ~isfield(options,'color')
        options.color = 'RGB';
    end
    
    if ~isfield(options,'channels')
        channels = [1 2 3]; % 1 is the first non-DAPI
    else
        channels = options.channels - 1; % in options.channels 1 = DAPI
    end
    if ~isfield(options,'FontSize')
        fs = 46;
    else
        fs = options.FontSize;
    end
    if ~isfield(options,'flip')
        options.flip = false;
    end

    pos = positions;
    radpos2Dcomb = radialPositive2D(pos, meta, channelThresholds, options);

    % now visualize
    bg = radpos2Dcomb.negative_norm_scaled;
    if ~isnan(channels(1))
        R = radpos2Dcomb.positive_norm_scaled{channels(1)};
    else
        R = bg*0;
    end
    if ~isnan(channels(2))
        G = radpos2Dcomb.positive_norm_scaled{channels(2)};
    else
        G = bg*0;
    end
    if ~isnan(channels(3))
        B = radpos2Dcomb.positive_norm_scaled{channels(3)};
    else
        B = bg*0;
    end

    datatable = [];
    datavars = {};
    
    total = radpos2Dcomb.anypositive_norm_scaled + bg;
    %total = radpos2Dcomb.density_scaled;
    %total = total/max(total(:));
    if strcmp(options.color, 'RGB')
        RGB = cat(3,R,G,B);
    elseif strcmp(options.color, 'CMY')
        RGB = cat(3,G+B, R+B, R+G);
	elseif strcmp(options.color, 'MG')
        RGB = cat(3,R, G, R);
    else
        error('invalid color mode');
    end
    RGBwmask = (max(total(:)) - total + RGB)/max(total(:));

    h = size(RGB,1);
    w = size(RGB,2);
    graphh = 260;
    pmarg = 130;
    lw = 7;
    d = 20; % additional offset for upper right
    figwidth = w + 2*pmarg + d + graphh;
    figheight = graphh+2*pmarg+h;
    
    figure('Position',[0 0 figwidth, figheight]),
    
    ax0 = subplot(2,2,1);

    Rmax = radpos2Dcomb.Rmax;
    Zmax = radpos2Dcomb.Zmax;
    if options.half
        imshow(RGBwmask,[],'Xdata',[0,Rmax],'Ydata',[0 Zmax])
    else
        imshow(RGBwmask,[],'Xdata',[-Rmax,Rmax],'Ydata',[0 Zmax])
    end

    % plot density outlines
    %-------------------------------
    hold on
    
    % individual colonies 
    colcol = lines(4);
    radpos2D = {}; 
    
    datamatrix = [];
    
    for pi = 1:numel(pos)

        options.Rmax = radpos2Dcomb.Rmax;
        options.Zmax = radpos2Dcomb.Zmax;
        radpos2D{pi} = radialPositive2D(pos(pi), meta, channelThresholds, options);

        % density outline for individual colonies    
        B = radpos2D{pi}.boundaries{1};
        if ~options.half
            B(:,2) = B(:,2) - Rmax;
        else
            B = B(B(:,1) > 1 & B(:,2) > 1,:);
        end
        plot(ax0, B(:,2), B(:,1),'-','LineWidth',2,'Color',[1 1 1]*0.5);%colcol(pi,:));

        
        % make a datamatrix to save everything to a single csv at the end
        datavars = [datavars, ['Col ' num2str(pi) ' bdry (x)'], ['Col ' num2str(pi) ' bdry (y)']];

        beforeRows = size(datamatrix,1);
        
        datamatrix(1:size(B,1),2*(pi-1)+1) = B(:,1)';
        datamatrix(1:size(B,1),2*pi) = B(:,2)';
        datamatrix(size(B,1)+1:end,2*(pi-1)+1:2*pi) = NaN; % pad these cols if shorter than previous table
        datamatrix(beforeRows+1:end,1:2*(pi-1)) = NaN; % pad other cols if this one is longer than previous table
    end
    
    % plot overall density outline
    B = radpos2Dcomb.boundaries{1};
    Rcutoff = max(B(:,2))*meta.xres;
    Zcutoff = max(B(:,1))*meta.xres;
    if ~options.half
        B(:,2) = B(:,2) - Rmax;
    else
        B = B(B(:,1) > 1 & B(:,2) > 1,:);
    end
    plot(ax0, B(:,2), B(:,1),'-','LineWidth',lw,'Color',[1 1 1]*0.5)
    %plot(B(:,2), B(:,1),'-','LineWidth',3,'Color','k')
    set(gca,'YDir','normal')
    
    datavars = [datavars, 'Combined bdry (x)', 'Combined bdry (y)'];

    beforeRows = size(datamatrix,1);
    beforeCols = size(datamatrix,2);
    datamatrix(1:size(B,1),beforeCols+1) = B(:,1)';
    datamatrix(1:size(B,1),beforeCols+2) = B(:,2)';
    datamatrix(size(B,1)+1:end,end-2:end) = NaN;
    datamatrix(beforeRows+1:end,1:beforeCols) = NaN;

    % % scatter plot some points
    % Npts = 5000;
    % R = radpos2Dcomb.R;
    % Z = radpos2Dcomb.Z;
    % idx = 1+round(rand([1 Npts])*(numel(R)-1));
    % %idx = 1:numel(R);
    % scatter(R(idx), Z(idx),1,'.k')

    %mypi = 1:numel(pos);
    % scattercolors = lines(4);
    % for pi = mypi
    % 
    %     Zp = pos(pi).cellData.Z*scalefactor;
    %     Rp = sqrt(sum((pos(pi).cellData.XY - pos(pi).center).^2,2));    
    %     scatter(Rp,Zp,'.','MarkerEdgeColor',scattercolors(pi,:));
    % end

    hold off

    % plot marginals: profiles in r
    %-------------------------------
    ax1 = subplot(2,2,3);
    if strcmp(options.color, 'RGB')
        %colors =[0.8 0 0; 0 0.8 0; [0 0 0.8]; 0 0 0];
        colors = lines(5);
        colors = cat(1, colors([2 5 1],:), [0 0 0]);
    elseif strcmp(options.color, 'CMY')
        colors =[0 1 1; 1 0 1; [1 0.8 0]; 0 0 0];
    elseif strcmp(options.color, 'MG')
        colors = lines(5);
        colors =[1 0 1; colors(5,:); 0 0 0; 0 0 0];
    end
    summed_profiles = {};
    for ci = 1:4
        summed_profiles{ci} = [];
    end

    hold on 

    r = radpos2D{1}.r_scaled*meta.xres;
    for pi = 1:numel(pos)

        for ci = 1:numel(radpos2D{pi}.positiveR_norm_scaled)
            sumprof = radpos2D{pi}.positiveR_norm_scaled{ci};
            sumprof(r > Rcutoff) = NaN;
            summed_profiles{ci} = cat(1, summed_profiles{ci}, sumprof);
            %plot(r, sumprof,'Color',colors(ci,:))
        end
        sumprof = radpos2D{pi}.negativeR_norm_scaled;
        sumprof(r > Rcutoff) = NaN;
        summed_profiles{4} = cat(1, summed_profiles{4}, sumprof);
        %plot(r, sumprof,'k');
    end

    summean = {};
    sumstd = {};
    channels = [channels 4];
    legendstr = {};
    subset = [];
    chanlabels = [meta.channelLabel, 'negative'];
    
    datavars = [datavars, 'r'];
    datamatrix(1:size(r',1),size(datamatrix,2)+1) = r';
    datamatrix(size(r',1)+1:end,end-1:end) = NaN;
    
    for cii = 1:4
        
        ci = channels(cii);
        
        if ~isnan(ci)
            
            summean{ci} = nanmean(summed_profiles{ci});
            sumstd{ci} = nanstd(summed_profiles{ci});
            
            good = ~isnan(summean{ci} + sumstd{ci});
            rdouble = [r(good), fliplr(r(good))];
            errdouble = [summean{ci}(good) + sumstd{ci}(good), fliplr(summean{ci}(good) - sumstd{ci}(good))];

            fill(rdouble, errdouble, colors(cii,:),'FaceAlpha',0.2,'EdgeColor','none');
            g = plot(r, summean{ci},'Color',colors(cii,:),'LineWidth',lw);
            
            legendstr = [legendstr, chanlabels(ci+1)];
            subset = [subset, g];

            % add to datamatrix to save everything to a single csv at the end
            datavars = [datavars, [chanlabels{ci+1} ' (r mean)'], [chanlabels{ci+1} ' (r std)']];

            beforeRows = size(datamatrix,1);
            beforeCols = size(datamatrix,2);
            
            datamatrix(1:size(summean{ci}',1),beforeCols+1) = summean{ci}';
            datamatrix(1:size(sumstd{ci}',1),beforeCols+2) = sumstd{ci}';
            datamatrix(size(summean{ci}',1)+1:end,end-2:end) = NaN;

            datamatrix(beforeRows+1:end,1:end-2) = NaN;
        end
    end

    hold off
    ylim([0 1])
    ylabel('positive fraction');
    %ylabel('P(+ | r)');
    xlabel('r (\mum)');

    if options.half
        xlim([0 radpos2Dcomb.Rmax]*meta.xres)
    else
        xlim([-radpos2Dcomb.Rmax radpos2Dcomb.Rmax]*meta.xres)
    end
    %legend(subset,legendstr,'Orientation','horizontal','Position',[0.35 0.4 0.1 0.2]);
    legend(subset,legendstr,'Orientation','horizontal','Position',[0.35 0.85 0.1 0.2]);
    legend('boxoff')
    
%     if options.half
%         xlim([0 Rcutoff + 10])
%     else
%         xlim([-(Rcutoff + 10) Rcutoff + 10])
%     end
    cleanSubplot(fs, lw)


    % plot marginals: profiles in z
    %-------------------------------

    % for i = 1:numel(channels)
    %         %errorbar(xi*meta.xres,P(channels(i),:),Pstd(channels(i),:),'LineWidth',lw,'Color',colors(i,:));
    %         good = ~isnan(P(channels(i),:));
    %         fill([xi,fliplr(xi)],...
    %             [P(channels(i),good) + Pstd(channels(i),good), fliplr(P(channels(i),good) - Pstd(channels(i),good))],...
    %             colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    % end

    ax2 = subplot(2,2,2);

    hold on 

    summed_profiles = {};
    for ci = 1:4
        summed_profiles{ci} = [];
    end

    z = radpos2D{1}.z_scaled'*meta.xres;
    for pi = 1:numel(pos)

        for ci = 1:3
            sumprof = radpos2D{pi}.positiveZ_norm_scaled{ci}';
            sumprof(z > Zcutoff) = NaN;
            summed_profiles{ci} = cat(1, summed_profiles{ci}, sumprof);
            %plot(r, sumprof,'Color',colors(ci,:))
        end
        sumprof = radpos2D{pi}.negativeZ_norm_scaled';
        sumprof(z > Zcutoff) = NaN;
        summed_profiles{4} = cat(1, summed_profiles{4}, sumprof);
    end

    summean = {};
    sumstd = {};
    
    datavars = [datavars, 'z'];
    datamatrix(1:size(z',1),size(datamatrix,2)+1) = z';
    datamatrix(size(z',1)+1:end,end-1:end) = NaN;
    
    for cii = 1:4
        
        ci = channels(cii);
        
        if ~isnan(ci)
            summean{ci} = nanmean(summed_profiles{ci});
            sumstd{ci} = nanstd(summed_profiles{ci});

            good = ~isnan(summean{ci} + sumstd{ci});
            zdouble = [z(good), fliplr(z(good))];
            errdouble = [summean{ci}(good) + sumstd{ci}(good), fliplr(summean{ci}(good) - sumstd{ci}(good))];

            fill(errdouble, zdouble, colors(cii,:),'FaceAlpha',0.2,'EdgeColor','none');
            plot(summean{ci}, z,'Color',colors(cii,:),'LineWidth',lw);
            
            % add to datamatrix to save everything to a single csv at the end
            datavars = [datavars, [chanlabels{ci+1} ' (z mean)'], [chanlabels{ci+1} ' (z std)']];

            beforeRows = size(datamatrix,1);
            beforeCols = size(datamatrix,2);
            
            datamatrix(1:size(summean{ci}',1),beforeCols+1) = summean{ci}';
            datamatrix(1:size(sumstd{ci}',1),beforeCols+2) = sumstd{ci}';
            datamatrix(size(summean{ci}',1)+1:end,end-2:end) = NaN;

            datamatrix(beforeRows+1:end,1:end-2) = NaN;
        end
    end

    hold off
    xlim([0 1])
    ylim([0 0.9*Zmax*meta.xres]);
    xlabel('positive fraction');
    %xlabel('P(+ | z)');
    ylabel('z (\mum)');
    
    cleanSubplot(fs,lw)
    
    % return data underlying graphs
    datatable = array2table(datamatrix, 'VariableNames',datavars);
    radpos2Dcomb.dataTable = datatable;
    
    % overall layout
    %---------------

    set(ax1 ,'Layer', 'Top')
    set(ax2 ,'Layer', 'Top')
    
    %linkaxes([ax0,ax1],'x');
    %linkaxes([ax0,ax2],'y');

    ax0.Units = "pixels";
    ax1.Units = "pixels";
    ax2.Units = "pixels";

    % Position: left bottom width height
    if options.flip
        ax0.Position = [pmarg pmarg w h]; % upper left
        ax1.Position = [pmarg h+3*pmarg/2 w graphh]; % lower left (origin)
        ax2.Position = [w+3*pmarg/2+d pmarg graphh h]; % upper right
    else
        ax0.Position = [pmarg graphh+3*pmarg/2 w h]; % upper left
        ax1.Position = [pmarg pmarg w graphh]; % lower left (origin)
        ax2.Position = [w+3*pmarg/2+d graphh+3*pmarg/2 graphh h]; % upper right
    end
end