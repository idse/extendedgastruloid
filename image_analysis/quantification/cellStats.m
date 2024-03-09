classdef cellStats < handle
    % Class to compare cell fate statistics

    % ---------------------
    % Idse Heemskerk, 2017
    % ---------------------

    properties

        % meta
        area
        conditions
        channelLabel
        markerChannels
        normalizeChannel

        % meta for Colony objects
        radiusPixel
        radiusMicron

        % basic stats
        nucLevel
        cytLevel
        XY
        sample
        lim
        scale               % scale = (max(A) - min(A));
        offset              % offset = min(A);
        bins
        histograms          % Nfiles x Nchannels
        othercelldata

        % clusters
        clustermodel
        clusters
        thresholds
        nClusters
        clusterHistograms
        cutoff
    end

    methods

        % constructor
        function this = cellStats(allPositions, meta, markerChannels, normalizeChannel)

            if exist('normalizeChannel','var')
                this.normalizeChannel = normalizeChannel;
            else
                this.normalizeChannel = [];
            end

            this.conditions = meta.conditions;
            this.channelLabel = meta.channelLabel;
            this.markerChannels = markerChannels;

            % area in cm^2
            this.area = meta.xSize*meta.ySize*meta.xres*meta.yres/10^8;

            % input is cell array
            if iscell(allPositions)

                if numel(allPositions) ~= meta.conditions
                    error('allData elements dont match number of conditions');
                end

                % make table of nuclear values
                for i = 1:numel(allPositions)

                    for ti = 1:meta.nTime
                        this.nucLevel{i,ti} = allPositions{i}.cellData(ti).nucLevel - allPositions{i}.cellData(ti).background;
                        if isfield(allPositions{i}.cellData(ti),'cytLevel')
                            this.cytLevel{i,ti} = allPositions{i}.cellData(ti).cytLevel - allPositions{i}.cellData(ti).background;
                        end
                    end
                end

                % input is Colony array
            elseif isa(allPositions,'Colony')%strcmp(class(allPositions),'Colony')

                posPerCondition = meta.posPerCondition;
                positions = allPositions;

                % combine data from the same condition
                for wellnr = 1:meta.nWells

                    if meta.loop4well
                        flipwell = 9-wellnr;
                        conditionPositions = [(posPerCondition/2)*(wellnr-1)+1:(posPerCondition/2)*wellnr...
                            (posPerCondition/2)*(flipwell-1)+1:(posPerCondition/2)*flipwell];
                        warning('loop4well will not work for unequal number of positions per condition');
                    else
                        conditionPositions = meta.conditionStartPos(wellnr):meta.conditionStartPos(wellnr)+meta.posPerCondition(wellnr)-1;
                        %conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
                    end

                    this.radiusPixel{wellnr} = positions(conditionPositions(1)).radiusPixel;
                    this.radiusMicron{wellnr} = positions(conditionPositions(1)).radiusMicron;

                    for ti = 1:meta.nTime
                        nucLevelCombined = [];
                        XYCombined = [];
                        sampleIdx = [];
                        for i = conditionPositions
                            nucLevel = positions(i).cellData(ti).nucLevel - positions(i).cellData(ti).background;
                            XYrel = positions(i).cellData(ti).XY - positions(i).center(ti,:);
                            XYrel = [XYrel, positions(i).cellData(ti).Z];
                            nucLevelCombined = cat(1,nucLevelCombined,nucLevel);
                            XYCombined = cat(1,XYCombined, XYrel);
                            sampleIdx = cat(1, sampleIdx, zeros([positions(i).ncells(ti) 1]) + i);
                        end
                        this.nucLevel{wellnr,ti} = nucLevelCombined;
                        this.XY{wellnr,ti} = XYCombined;
                        this.sample{wellnr, ti} = sampleIdx;

                        if isfield(positions(1).cellData,'cytLevel')
                            cytLevelCombined = [];
                            for i = conditionPositions
                                cytLevel = positions(i).cellData(ti).cytLevel - positions(i).cellData(ti).background;
                                cytLevelCombined = cat(1,cytLevelCombined,cytLevel);
                            end
                            this.cytLevel{wellnr,ti} = cytLevelCombined;
                        end
                    end

                    % all the other properties
                    propnames = setdiff(fields(positions(1).cellData),{'cytLevel','nucLevel','XY'});
                    for ti = 1:meta.nTime
                        this.othercelldata{wellnr,ti} = struct();
                        for propi = 1:numel(propnames)
                            propvalcombined = [];
                            for i = conditionPositions
                                propval = positions(i).cellData(ti).(propnames{propi});
                                propvalcombined = cat(1,propvalcombined,propval);
                            end
                            this.othercelldata{wellnr,ti}.(propnames{propi}) = propvalcombined;
                        end
                    end
                end

            elseif isa(allPositions,'Position')%strcmp(class(allPositions),'Position')

                posPerCondition = meta.posPerCondition;
                positions = allPositions;

                % combine data from the same condition
                for wellnr = 1:meta.nWells

                    if meta.loop4well
                        flipwell = 9-wellnr;
                        conditionPositions = [(posPerCondition/2)*(wellnr-1)+1:(posPerCondition/2)*wellnr...
                            (posPerCondition/2)*(flipwell-1)+1:(posPerCondition/2)*flipwell];
                        warning('loop4well will not work for unequal number of positions per condition');
                    else
                        conditionPositions = meta.conditionStartPos(wellnr):meta.conditionStartPos(wellnr)+meta.posPerCondition(wellnr)-1;
                        %conditionPositions = posPerCondition*(wellnr-1)+1:posPerCondition*wellnr;
                    end
                    
                    for ti = 1:meta.nTime
                        nucLevelCombined = [];
                        XYCombined = [];
                        sampleIdx = [];
                        for i = conditionPositions
                            nucLevel = positions(i).cellData(ti).nucLevel - positions(i).cellData(ti).background;
                            %                             XYrel = positions(i).cellData(ti).XY - positions(i).center(ti,:);
                            nucLevelCombined = cat(1,nucLevelCombined,nucLevel);
                            %                             XYCombined = cat(1,XYCombined, XYrel);
                            sampleIdx = cat(1, sampleIdx, zeros([positions(i).ncells(ti) 1]) + i);
                        end
                        this.nucLevel{wellnr,ti} = nucLevelCombined;
                        this.XY{wellnr,ti} = XYCombined;
                        this.sample{wellnr, ti} = sampleIdx;

                        if isfield(positions(1).cellData,'cytLevel')
                            cytLevelCombined = [];
                            for i = conditionPositions
                                cytLevel = positions(i).cellData(ti).cytLevel - positions(i).cellData(ti).background;
                                cytLevelCombined = cat(1,cytLevelCombined,cytLevel);
                            end
                            this.cytLevel{wellnr,ti} = cytLevelCombined;
                        end
                    end

                    % all the other properties
                    propnames = setdiff(fields(positions(1).cellData),{'cytLevel','nucLevel','XY'});
                    for ti = 1:meta.nTime
                        this.othercelldata{wellnr,ti} = struct();
                        for propi = 1:numel(propnames)
                            propvalcombined = [];
                            for i = conditionPositions
                                propval = positions(i).cellData(ti).(propnames{propi});
                                propvalcombined = cat(1,propvalcombined,propval);
                            end
                            this.othercelldata{wellnr,ti}.(propnames{propi}) = propvalcombined;
                        end
                    end
                end

            end

            % optional normalization by some channel (think DAPI)
            for i = 1:numel(this.nucLevel)
                if ~isempty(this.normalizeChannel)
                    N = this.nucLevel{i}(:,normalizeChannel);
                    this.nucLevel{i} = bsxfun(@rdivide, this.nucLevel{i}, N);
                end
            end
        end

        function setLimits(this, tol)
            % determine intensity limits for distributions and plots

            nucLevelAll = [];
            cytLevelAll = [];
            for i = 1:numel(this.nucLevel)
                nucLevelAll = cat(1, nucLevelAll, this.nucLevel{i});
                if ~isempty(this.cytLevel)
                    cytLevelAll = cat(1, cytLevelAll, this.cytLevel{i});
                end
            end

            this.lim = {};
            this.scale = zeros([1 4]);
            this.offset = zeros([1 4]);
            for i = 1:size(nucLevelAll, 2)
                if ~isempty(cytLevelAll)
                    A = cat(1,nucLevelAll(:,i),cytLevelAll(:,i));
                else
                    A = nucLevelAll(:,i);
                end
                this.scale(i) = (max(A) - min(A));
                this.offset(i) = min(A);
                this.lim{i} = stretchlim(mat2gray(A),tol)*(max(A) - min(A)) + min(A);
                this.lim{i}(1) = max(min(A),0);
                if min(A) < 0
                    str = [num2str(sum(A < 0)) ' out of ' num2str(numel(A))];
                    warning([str ' values below background found']);
                end
            end

            %             for i = 1:numel(this.nucLevel)
            %                 this.lim{i} = this.lim{i} - this.offset(i);
            %                 for j = 1:size(nucLevelAll, 2)
            %                     this.nucLevel{i}(:,j) = this.nucLevel{i}(:,j) - this.offset(j);
            %                 end
            %             end
        end

        function makeHistograms(this, nBins, tolerance)
            % make histograms
            %
            % makeHistograms(nBins)

            % set intensity limits
            this.setLimits(tolerance);

            this.histograms = {};
            this.bins = {};

            nConditions = size(this.nucLevel, 1);
            nTime = size(this.nucLevel, 2);
            for ti = 1:nTime
                for fi = 1:nConditions
                    for channelIndex = 1:numel(this.lim)

                        this.bins{channelIndex} = linspace(this.lim{channelIndex}(1), this.lim{channelIndex}(2), nBins);
                        nn = histc(this.nucLevel{fi,ti}(:,channelIndex), this.bins{channelIndex});
                        this.histograms{fi, channelIndex, 1, ti} = nn;

                        if ~isempty(this.cytLevel) && ~isempty(this.cytLevel{fi})
                            nc = histc(this.cytLevel{fi,ti}(:,channelIndex), this.bins{channelIndex});
                            this.histograms{fi, channelIndex, 2, ti} = nc;
                        end
                    end
                end
            end
        end

        function exportCSV(this, meta)

            for condi = 1:numel(meta.conditions)
                disp(meta.conditions{condi})
                M = [this.nucLevel{condi}, this.XY{condi}, this.sample{condi}];
                T = array2table(M);
                T.Properties.VariableNames =  [meta.channelLabel {'X','Y','sample'}];
                fname = [meta.conditions{condi} '_' [meta.channelLabel{2:end}] '.csv'];
                writetable(T,fullfile(dataDir, fname))
            end
        end

        function thresh = getThresholds(this, confidence, conditions, whichthreshold)

            if ~exist('confidence','var')
                confidence = 0.95;
            end
            if ~exist('conditions','var')
                conditions = 1:numel(this.nucLevel);
            end
            if ~exist('whichthreshold','var')
                whichthreshold = [];%ones([1 max(this.markerChannels)]);
            end
            if ~exist('SharedCovariance','var')
                SharedCovariance = 0;
            end

            disp('determining thresholds')
            tic
            % fit Gaussian mixture model to each channel separately
            BIC = zeros([6 max(this.markerChannels)]);

            kmax = 6;
            thresh = zeros([kmax max(this.markerChannels)]);
            regularizationVal = 0;

            %X = stats.nucLevel{condi};
            X = cat(1,this.nucLevel{conditions});
            disp('setting vals below background to background');
            X(X < 0) = 0;
            nucInorm = log(1 + X./mean(X,'omitnan'));

            for k = 1:kmax

                for ci = this.markerChannels

                    %y = (nucInorm(:,ci) > mean(nucInorm(:,ci)) + std(nucInorm(:,ci))) + 1;
                    %        'Start',y,...

                    try
                        clustermodel1d = fitgmdist(nucInorm(:,ci),k,'CovarianceType','full','SharedCovariance',false,...
                            'RegularizationValue', regularizationVal,'Replicates',1,...
                            'Options',statset('MaxIter',500,'TolFun',1e-6,'Display','off'));
                    catch
                        clustermodel1d = fitgmdist(nucInorm(:,ci),k,'CovarianceType','full','SharedCovariance',true,...
                            'RegularizationValue', regularizationVal,'Replicates',1,...
                            'Options',statset('MaxIter',500,'TolFun',1e-6,'Display','off'));
                        warning('Ill-conditioned covariance created, use SharedCovariance');
                    end

                    [clusterLabels, ~, clusterProb, ~] = cluster(clustermodel1d, nucInorm(:,ci));
                    BIC(k, ci) = clustermodel1d.BIC;

                    meanvec = zeros([k 1]);
                    for clusti = 1:k
                        meanvec(clusti) = mean(X(clusterLabels==clusti,ci));
                    end

                    % the negative cluster label (one with the smaller mean
                    [~,clusteridx] = sort(meanvec);
                    %[~,negcluster] = min(meanvec);
                    if isempty(whichthreshold) && k > 1 % by default only the highest mean cluster is positive
                        negcluster = clusteridx(end-1);
                    elseif ~isempty(whichthreshold) && whichthreshold(ci) <= k
                        negcluster = clusteridx(whichthreshold(ci));
                    else
                        negcluster = clusteridx(1);
                    end
                    idx = clusterProb(:,negcluster) > 1-confidence;

                    thresh(k,ci) = max(X(idx,ci));
                end
            end

            BICnorm = (BIC-min(BIC,[],1))./(max(BIC,[],1)-min(BIC,[],1));
            % this determines the smallest k for which the BIC is within 10% of its
            % minimum
            [~,kvals] = min(BICnorm.*(BICnorm > 0.1),[],1);

            disp(['k : ' num2str(kvals)])

            finalthresh = zeros([1 max(this.markerChannels)]);
            for ci = this.markerChannels
                finalthresh(ci) = thresh(kvals(ci),ci);
            end
            this.thresholds = finalthresh;
            toc

            % diagnosis
            %             ci = 2;
            %             cj = 3;
            %             k = 3;
            %             scatter(nucInorm(:,ci),nucInorm(:,cj),2,nucInorm(:,ci) > log(1+thresh(k,ci)/mean(X(:,ci))) )

            % visualize marginal distribution of one marker
            %'Bandwidth',bw
            %[fpos,xi] = ksdensity(nucInorm(:,i),'BoundaryCorrection','reflection','Support','positive');
            %plot(xi,fpos)

            % Otsu alternative (doesn't work well)
            %             condi = 1;
            %             for ci = 2:4
            %                 Y = stats.nucLevel{condi}(:,ci);
            %                 X = log(1 + Y/mean(Y));
            %                 t = otsuthresh(X);
            %                 A = max(X)-min(X);
            %                 ot(ci) = mean(Y)*(exp(min(X) + t*A)-1);
            %             end
        end

        %-----------------------------------------------------------------
        % clusters OLD NEEDS REVISIONS
        %-----------------------------------------------------------------

        function makeClustersManual(this, k, ci1, ci2)

            this.nClusters = k;
            this.cutoff = 0.5;

            % combine nuclear values of all conditions
            nucI = this.nucLevel{1};
            condIdx = ones([size(nucI,1) 1],'uint8');
            for fi = 2:numel(this.conditions)
                nucI = cat(1, nucI, this.nucLevel{fi});
                condIdx = cat(1, condIdx, fi*ones([size(this.nucLevel{fi},1) 1],'uint8'));
            end

            figure,
            scatter(nucI(:,ci1), nucI(:,ci2),1)
            xlim(this.lim{ci1}');
            ylim(this.lim{ci2}');
            xlabel(this.channelLabel{this.markerChannels(ci1)});
            ylabel(this.channelLabel{this.markerChannels(ci2)});


            clusterLabels = 0*nucI(:,ci1);
            clusterProb = 0*nucI(:,ci1);
            for i = 1:k
                h(i) = impoly;
                clusterX = h(i).getPosition;
                idx = inpolygon(nucI(:,ci1),nucI(:,ci2),clusterX(:,1),clusterX(:,2));
                clusterLabels(idx) = i;
                clusterProb(idx) = 1;
                this.clusters{i} = struct();
            end
            close;

            % cluster every individual conditions
            this.clusters = {};
            for fi = 1:numel(this.conditions)

                idx = condIdx == fi;
                this.clusters{fi} = struct('clusterLabels', clusterLabels(idx),...
                    'clusterProb',  clusterProb(idx,:));

                % add mean and std
                for k = 1:this.nClusters
                    clusterIdx = this.getClusterIdx(fi, k);
                    this.clusters{fi}.mean{k} = mean(this.nucLevel{fi}(clusterIdx,:));
                    this.clusters{fi}.std{k} = std(this.nucLevel{fi}(clusterIdx,:));
                end
            end

        end

        function [clusterLabels, nucInorm] = makeClusters(this, k, regularizationVal, cutoff)

            this.cutoff = cutoff;
            this.nClusters = k;

            mc = this.markerChannels;

            % combine nuclear values of all conditions
            nucI = this.nucLevel{1}(:,mc);
            condIdx = ones([size(nucI,1) 1],'uint8');
            for fi = 2:numel(this.conditions)
                nucI = cat(1, nucI, this.nucLevel{fi}(:,mc));
                condIdx = cat(1, condIdx, fi*ones([size(this.nucLevel{fi}(:,mc),1) 1],'uint8'));
            end

            % remove outliers and normalize the values
            nucInorm = nucI;
            for i = 1:numel(mc)
                outliers = nucI(:,i) > this.lim{mc(i)}(2);
                %nucI(outliers,:) = [];
                %nucInorm(outliers,:) = [];
                nucInorm(:,i) = nucInorm(:,i) - this.lim{mc(i)}(1);
                nucInorm(:,i) = nucInorm(:,i)./max(nucInorm(~outliers,i));
            end

            % fit Gaussian mixture model to all data
            this.clustermodel = fitgmdist(nucInorm(~outliers,:),k,'CovarianceType','full','SharedCovariance',false,...
                'RegularizationValue', regularizationVal);

            % cluster all
            [clusterLabels, ~, clusterProb, ~] = cluster(this.clustermodel, nucInorm);
            clusterLabels(outliers) = 0;

            % cluster every individual conditions
            this.clusters = {};
            for fi = 1:numel(this.conditions)

                idx = condIdx == fi;
                this.clusters{fi} = struct('clusterLabels', clusterLabels(idx),...
                    'clusterProb',  clusterProb(idx,:));

                % add mean and std
                for k = 1:this.nClusters
                    clusterIdx = this.getClusterIdx(fi, k);
                    this.clusters{fi}.mean{k} = mean(this.nucLevel{fi}(clusterIdx,:));
                    this.clusters{fi}.std{k} = std(this.nucLevel{fi}(clusterIdx,:));
                end
            end

            %             % cluster every individual conditions, by applying model to
            %             each case separately -> didn't work that well
            %             this.clusters = {};
            %             for fi = 1:numel(this.conditions)
            %
            %                 Xc = this.nucLevel{fi}(:,this.markerChannels);
            %                 Xcnorm = Xc;
            %                 for i = 1:numel(this.markerChannels)
            %                     Xcnorm(:,i) = Xcnorm(:,i) - this.lim{this.markerChannels(i)}(1);
            %                     Xcnorm(:,i) = Xcnorm(:,i)./max(Xcnorm(:,i));
            %                 end
            %
            %                 [clusterLabels, ~, clusterProb, ~] = cluster(this.clustermodel, Xcnorm);
            %                 this.clusters{fi} = struct('clusterLabels', clusterLabels, 'clusterProb',  clusterProb);
            %
            %                 % add mean and std
            %                 for k = 1:this.nClusters
            %                     clusterIdx = this.getClusterIdx(fi, k);
            %                     this.clusters{fi}.mean{k} = mean(this.nucLevel{fi}(clusterIdx,:));
            %                     this.clusters{fi}.std{k} = std(this.nucLevel{fi}(clusterIdx,:));
            %                 end
            %             end
        end

        function clusterIdx = getClusterIdx(this, conditionIdx, clusterLabel)

            C = this.clusters{conditionIdx};
            clusterIdx = C.clusterLabels == clusterLabel;
            goodIdx = max(C.clusterProb,[],2) > this.cutoff;
            clusterIdx = clusterIdx & goodIdx;
        end

        function makeClusterHistograms(this)

            for clusterLabel = 1:this.nClusters
                for fi = 1:numel(this.conditions)

                    clusterIdx = this.getClusterIdx(fi, clusterLabel);

                    for channelIndex = 1:numel(this.lim)

                        n = histc(this.nucLevel{fi}(clusterIdx,channelIndex), this.bins{channelIndex});
                        if sum(n) < 10
                            n = zeros([numel(this.bins{channelIndex}) 1],class(n));
                        end
                        this.clusterHistograms{fi, channelIndex, clusterLabel} = n;
                    end
                end
            end
        end

        %-----------------------------------------------------------------
        % summary
        %-----------------------------------------------------------------

        function printSummary(this)

            clusterLabel = [];
            intensitySummary(clusterLabel);

            %-----------------------------------

            if ~isempty(this.clusters)

                for clusterLabel = 1:this.nClusters
                    intensitySummary(clusterLabel)
                end

                %-----------------------------------

                width = 12;
                conditionstr = this.conditions;
                for i = 1:numel(this.conditions)
                    N = numel(conditionstr{i}) + 6;
                    conditionstr{i} = strpad(conditionstr{i}, N,'post',' ');
                end
                t = [strpad('', width,'post',' ') conditionstr{:}];

                disp(' ');
                disp(t)
                disp(repmat('-',1,numel(t)));

                % cell numbers
                s = strpad(['Ncells '], width,'post',' ');
                valstr = [];
                for fi = 1:numel(this.conditions)

                    N = numel(conditionstr{fi});
                    Nt = size(this.nucLevel{fi},1);
                    valstr = [valstr strpad(sprintf('%d', Nt), N,'post',' ')];
                end
                disp([s valstr]);

                % density
                s = strpad('cells/cm^2', width,'post',' ');
                valstr = [];
                for fi = 1:numel(this.conditions)

                    N = numel(conditionstr{fi});
                    rho = size(this.nucLevel{fi},1)/this.area;
                    valstr = [valstr strpad(sprintf('%.1e', rho), N,'post',' ')];
                end
                disp([s valstr]);

                % fractions in cluster
                for clusterLabel = 1:this.nClusters

                    s = strpad(['cluster ' num2str(clusterLabel)], width,'post',' ');
                    valstr = [];
                    for fi = 1:numel(this.conditions)

                        N = numel(conditionstr{fi});
                        Nt = size(this.clusters{fi}.clusterLabels,1);
                        Nc = sum(this.getClusterIdx(fi,clusterLabel));
                        fraction = Nc/Nt;
                        valstr = [valstr strpad(sprintf('%.2f', fraction), N,'post',' ')];
                    end
                    disp([s valstr]);
                end
            end

            function intensitySummary(clusterLabel)

                width = 16;
                N = max(cellfun(@numel, this.conditions)) + 5;
                labels = [];
                for i = 1:numel(this.markerChannels)
                    labels = [labels strpad(this.channelLabel{this.markerChannels(i)}, width,'post',' ')];
                end
                if isempty(clusterLabel)
                    t = [strpad('',N,'post',' ') labels];
                else
                    t = [strpad(['cluster ' num2str(clusterLabel)], N,'post',' ') labels];
                end

                disp(' ');
                disp(t)
                disp(repmat('-',1,numel(t)));
                for fj = 1:numel(this.conditions)

                    if isempty(clusterLabel)
                        idx = true([size(this.nucLevel{fj},1) 1]);
                    else
                        idx = this.getClusterIdx(fj, clusterLabel);
                    end

                    if sum(idx) < 10
                        idx = idx*0 > 0;
                    end

                    s = strpad(this.conditions{fj}, N,'post',' ');

                    meanval = mean(this.nucLevel{fj}(idx, this.markerChannels));
                    stdval = std(this.nucLevel{fj}(idx, this.markerChannels));

                    valstr = [];
                    for channel = 1:numel(this.markerChannels)
                        valstr = [valstr strpad(sprintf('%.0f (%.0f)', meanval(channel), stdval(channel)), width,'post',' ')];
                    end
                    disp([s valstr])
                end
            end
        end

        %-----------------------------------------------------------------
        % visualization
        %-----------------------------------------------------------------

        function makeScatterPlot(this, conditionIdx, channelIdx, showClusters, cytoplasmic)

            if ~exist('showClusters','var')
                showClusters = false;
            end
            if ~exist('cytoplasmic','var')
                cytoplasmic = [false false];
            end

            s = 50;
            c1 = channelIdx(1);
            c2 = channelIdx(2);

            cdi = conditionIdx;

            hold on

            if cytoplasmic(1)
                A = this.cytLevel{cdi}(:,c1);
                xlabelstr = [this.channelLabel{c1} ' cytoplasmic'];
            else
                A = this.nucLevel{cdi}(:,c1);
                xlabelstr = this.channelLabel{c1};
            end
            if cytoplasmic(2)
                B = this.cytLevel{cdi}(:,c2);
                ylabelstr = [this.channelLabel{c2} ' cytoplasmic'];
            else
                B = this.nucLevel{cdi}(:,c2);
                ylabelstr = this.channelLabel{c2};
            end

            colormap(lines(this.nClusters));
            %colors = lines(this.nClusters);
            fs = 20;

            if ~showClusters
                scatter(A, B, s,'.');
                disp('no clusters');
            else
                clusterLabels = this.clusters{cdi}.clusterLabels;
                idx = max(this.clusters{cdi}.clusterProb,[],2) > this.cutoff;
                scatter(A, B, s, '.','MarkerEdgeColor', 0.6*[1 1 1]);

                scatterlegend = {'unassigned'};

                for i = 1:max(clusterLabels)
                    idxx = clusterLabels == i & idx;
                    scatter(A(idxx), B(idxx), s, clusterLabels(idxx),'.'); %colors(clusterLabels(idxx),:)
                    scatterlegend = [scatterlegend, ['cluster ' num2str(i)]];
                end

                legend(scatterlegend)
            end

            xlabel(xlabelstr, 'FontSize',fs, 'FontWeight','Bold');
            ylabel(ylabelstr, 'FontSize',fs, 'FontWeight','Bold');
            title([this.conditions{cdi} ' (corr ' num2str(corr2(A,B),2) ')'],...
                'Interpreter','none', 'FontSize',fs, 'FontWeight','Bold');

            %axis([this.lim{c1}' this.lim{c2}']);

            set(gcf,'color','w');
            set(gca, 'LineWidth', 2);
            set(gca,'FontSize', fs)
            set(gca,'FontWeight', 'bold')

            hold off
        end

        %-----------------------------------------------------------------

        function plotDistributionComparison(this, options)

            if ~isfield(options,'channelIndex')
                error('please specify options.channelIndex');
            elseif numel(options.channelIndex) > 1
                error('please specify only one channel in options.channelIndex');
            else
                channelIndex = options.channelIndex;
            end
            if ~isfield(options,'cumulative')
                options.cumulative = false;
            end
            if ~isfield(options, 'clusterLabel')
                options.clusterLabel = [];
            end
            if ~isfield(options, 'conditionIdx')
                options.conditionIdx = 1:size(this.histograms,1);
            end
            if ~isfield(options,'time')
                ti = 1;
            else
                ti = options.time;
            end

            if isfield(options,'titlestr')
                titlestr = options.titlestr;
            elseif ~isempty(this.channelLabel)
                titlestr = this.channelLabel{channelIndex};
            else
                titlestr = ['channel ' num2str(channelIndex)];
            end

            if ~isfield(options,'cytoplasmic') || options.cytoplasmic == false
                j = 1;
            else
                j = 2;
                titlestr = [titlestr ' cytoplasmic'];
            end

            if ~isempty(options.clusterLabel)
                nall = [this.clusterHistograms{:, channelIndex, options.clusterLabel}];
                titlestr = [titlestr ';  cluster ' num2str(options.clusterLabel)];
            else
                nall = [this.histograms{options.conditionIdx, channelIndex, j, ti}];
            end

            if options.cumulative
                dist = bsxfun(@rdivide, nall, sum(nall,1));
                dist = cumsum(dist);
            else
                dist = bsxfun(@rdivide, nall, sum(nall,1));
                %dist = bsxfun(@rdivide, nall, max(nall,[],1));
            end

            nConditions = size(nall, 2);
            if nConditions < 8
                colors = lines(nConditions);
            else
                colors = hsv(nConditions);
            end

            [x,y] = histForBarlikePlot(this.bins{channelIndex}, dist);

            hold on
            fs = 20;
            for i = 1:nConditions
                plot(x,y(:,i),'LineWidth',2, 'Color',colors(i,:));
            end
            title(titlestr, 'FontSize',fs, 'FontWeight','Bold');
            hold off

            xlabel('Intensity', 'FontSize',fs, 'FontWeight','Bold');
            ylabel('Frequency', 'FontSize',fs, 'FontWeight','Bold');

            set(gcf,'color','w');
            set(gca, 'LineWidth', 2);
            set(gca,'FontSize', fs)
            set(gca,'FontWeight', 'bold')

            xlim(this.lim{channelIndex});

            if options.cumulative
                legend(this.conditions(options.conditionIdx), 'Location','SouthEast','FontSize',fs/2);
                ylim([0 1.05]);
            else
                legend(this.conditions(options.conditionIdx), 'Location','NorthEast','FontSize',fs/2);
                %ylim([0 0.3]);
                %ylim([0 1.05]);
            end

        end

        %-----------------------------------------------------------------

        function showClusters(this, channelIdx)
            % plot Gaussian covariance as ellipse
            % UNFINISHED
            % NEED TO STORE VALUES TO SCALE BACK TO ORIGINAL COORDINATES

            hold on
            fs = 15;
            gmfit = this.clustermodel;

            c1 = channelIdx(1);
            c2 = channelIdx(2);

            ci1 = this.markerChannels == c1;
            ci2 = this.markerChannels == c2;

            %         % means of gaussians
            %         scatter(gmfit.mu(:,1), gmfit.mu(:,2),100,'.k')

            % covariance ellipse and cluster labels
            for i = 1:this.nClusters

                text(gmfit.mu(i,ci1), gmfit.mu(i,ci2), num2str(i),'FontSize',fs, 'FontWeight','Bold');

                [V,D] = eig(gmfit.Sigma(:,:,i));
                t = linspace(0,2*pi);

                % major and minor axis length
                lmin = sqrt(D(1,1));
                lmaj = sqrt(D(2,2));

                majV = V(:,2);
                phi = atan(majV(2)/majV(1));
                R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

                %quiver(gmfit.mu(i,1), gmfit.mu(i,2), majV(1), majV(2));

                X = lmaj*cos(t);
                Y = lmin*sin(t);

                Xp = R(1,1)*X + R(2,1)*Y + gmfit.mu(i,ci1);
                Yp = R(1,2)*X + R(2,2)*Y + gmfit.mu(i,ci2);

                plot(Xp, Yp,'LineWidth',2,'Color','k')
            end

            hold off
        end
    end
end