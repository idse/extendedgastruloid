classdef Colony < Position
    % Data class to store a stem cell colony
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    properties

        center          % x-y pixel coordinates of center (relative to btf)
        well            % well located in
        
        % different from cellTracker: 
        %------------------------------
        
        radiusPixel     % radius of colony
        radiusMicron
        
        radialProfile   % struct with fields
                        % NucAvg    : radial average 
                        % NucAvgSeg : radial average based on segmentation
                        % NucStd    : radial standard deviation
                        % NucStdSeg  
                        % CytAvg    : cytoplasmic
                        % CytAvgSeg 
                        % CytStd
                        % CytStdSeg 
                        % BinEdges  : radial bin edges in pixels
                        % assumed to be the same for segmented and
                        % unsegmented
        
        %CM              % center of mass -> should use '.center'
        I               % second moments (inertia tensor)
    end
    
    properties (Dependent)
        radius;         % =radiusPixel for compatibility with CellTracker
    end
    
    methods
        
        % constructor
        function this = Colony(varargin)
            % Colony(meta)
            % Colony() (for initializing Colony arrays)
            
            % old: Colony(nChannels, center, radiusPixel, radiusMicron, boundingBox, well, nTime)

            this@Position(varargin{:});
            
            %this.center = center;
            %this.radiusPixel = radiusPixel;
            %this.radiusMicron = radiusMicron;
        end

        function debugInfo = extractData(this, dataDir, opts)
            
            debugInfo = extractData@Position(this, dataDir, opts);
            this.center = zeros([this.nTime 2]);
            
            this.setCenter();
        end
        
        function setCenter(this, margin)
            % setCenter()
            % setCenter(margin) % margin in microns
            
            if ~exist('margin','var')
                margin = 20;
            end
            
            for ti = 1:this.nTime
                
                % extractData setting center based on mean cell centroid
                CM = mean(this.cellData(ti).nucArea.*this.cellData(ti).XY)/mean(this.cellData(ti).nucArea);
                this.center(ti,:) = CM;
            
                % exclude cells/junk outside colony
                d = sqrt((this.cellData(ti).XY(:,1) - this.center(ti,1)).^2 ...
                        + (this.cellData(ti).XY(:,2) - this.center(ti,2)).^2);
                
                xyres = this.radiusMicron/this.radiusPixel;
                outside = d > this.radiusPixel + margin/xyres;

                fields = fieldnames(this.cellData(ti));
                fields = fields(cellfun(@(x) ~strcmp(x,'background'),fields));
                if ~isfield(this.cellData(ti),'cytLevel') || isempty(this.cellData(ti).cytLevel)
                    fields = fields(cellfun(@(x) ~strcmp(x,'cytLevel'),fields));
                end

                for i=1:numel(fields) % - 1 to exclude background which is not for each cell
                  this.cellData(ti).(fields{i}) = this.cellData(ti).(fields{i})(~outside,:);
                end
                
                this.ncells = sum(~outside);
                
                % recenter after removing junk outside
                CM = mean(this.cellData(ti).nucArea.*this.cellData(ti).XY)/mean(this.cellData(ti).nucArea);
                this.center(ti,:) = CM;
            end
        end
        
        
        % process data
        %---------------------------------

        function makeRadialAvgSeg(this, channels, normChannel, badidx)
            % create radial profile of segmented single cell data 
            % also calculate moment of inertia as measure of asymmetry
            % 
            % makeRadialAvgSeg()
            %
            % populates this.radialProfile.AvgSeg where rows correspond to
            % radial bins, and columns to channels
            
            warning('this function averages in radial bins that may contain very few points on the edge, leading to artefacts of averages going down on the edge, plotRadialProfiles no longer uses this as input but computes with fixed numbers of points per bin directly from stats');
            
            if ~exist('channels','var') || isempty(channels)
                channels = 1:numel(this.dataChannels);
                if isempty(this.dataChannels)
                    error('dataChannel and channels empty');
                end
            end
            
            this.radialProfile = struct('BinEdges',[],'NucAvgSeg',[],...
                'NucStdSeg',[],'CytAvgSeg',[],'CytStdSeg',[],'NucCytRatio',[],...
                'nCells',[],'NucArea',[],'CytMaskArea',[]);
            
            binWidthMicron = 8; % about two cell widths
            N = this.radiusMicron/binWidthMicron;
            binEdges = sqrt(linspace(0,this.radiusPixel^2,N+1));
            
            for ti = 1:numel(this.cellData)

                if ~exist('badidx','var') || isempty(badidx)
                    badidxti = false(size(this.cellData(ti).XY(:,1)));
                else
                    badidxti = badidx{ti};
                end
            
                if ~isempty(this.cellData(ti).XY)
    
                radialProf = struct('BinEdges',binEdges);
                
                XY = this.cellData(ti).XY;
                XY(:,1) = XY(:,1) - mean(XY(:,1));
                XY(:,2) = XY(:,2) - mean(XY(:,2));
                
                % determine which bin each cell is in
                r = sqrt(sum(XY(:,1:2).^2,2));
                [n,bini] = histc(r, radialProf.BinEdges);
                nBins = numel(n)-1; % last bin we ignore (see doc histc)
                
                N = numel(this.dataChannels);
                radialProf.NucAvgSeg = zeros([nBins N]);
                radialProf.NucStdSeg = zeros([nBins N]);
                radialProf.CytAvgSeg = zeros([nBins N]);
                radialProf.CytStdSeg = zeros([nBins N]);
                radialProf.nCells = zeros([nBins 1]);
                radialProf.NucArea = zeros([nBins 1]);
                radialProf.NucCytRatio = zeros([nBins N]);
                radialProf.CytMaskArea = zeros([nBins 1]);
                
                if isfield(this.cellData(ti),'cytLevel')
                    cyt = ~isempty(this.cellData(ti).cytLevel);
                else
                    cyt = false;
                end

                for cii = channels
                    
                    for i = 1:nBins
                        
                        % for the cells in bini calculate properties
                        idx = bini == i & ~badidxti;
                        nucbindata = this.cellData(ti).nucLevel(idx, cii)...
                                            - this.cellData(ti).background(cii);

                        % make nuc:cyt before potentially normalizing nuc
                        % by DAPI levels below
                        if cyt 
                            cytbindata = this.cellData(ti).cytLevel(idx, cii)...
                                                - this.cellData(ti).background(cii);
                            radialProf.CytAvgSeg(i,cii) = nanmean(cytbindata);
                            radialProf.CytStdSeg(i,cii) = nanstd(cytbindata);
                            
                            R = nucbindata./cytbindata;
                            %R(R < 0.3) = [];
                            %R(R > 2) = [];
                            radialProf.NucCytRatio(i,cii) = nanmean(R);
                        end
                        
                        % 180626 the behavior is different here from old code
                        % commented out below, where normChannel = 0
                        % normalized by cytoplasm, because here NucCytRatio
                        % is its own thing, this may affect Smad2 analysis
                        % if rerun
                        if exist('normChannel','var') && ~isempty(normChannel) && normChannel ~= 0
                            bg = this.cellData(ti).nucLevel(idx, normChannel)...
                                    - this.cellData(ti).background(normChannel);
                            nucbindata = nucbindata./bg;
                        end     
                        
                        radialProf.NucAvgSeg(i,cii) = nanmean(nucbindata);
                        radialProf.NucStdSeg(i,cii) = nanstd(nucbindata);
                        radialProf.nCells(i) = sum(idx);
                        radialProf.NucArea(i) = nanmean(this.cellData(ti).nucArea(idx));
                        if cyt
                            radialProf.CytMaskArea(i) = nanmean(this.cellData(ti).cytMaskArea(idx));
                        end
                    end
                end
                this.radialProfile(ti) = radialProf;
            end
            end
            
            % make radial time traces
            for cii = channels
                
                this.timeTraces.radAvgNuc{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgCyt{cii} = zeros([this.nTime nBins]);
                this.timeTraces.radAvgRatio{cii} = zeros([this.nTime nBins]);

                this.timeTraces.moment{cii} = zeros([this.nTime 2]);
                
                for ti = 1:numel(this.cellData)
                    
                    % moment
                    XY = this.cellData(ti).XY - mean(this.cellData(ti).XY);
                    nucdata = this.cellData(ti).nucLevel(:,cii);
                    this.timeTraces.moment{cii}(ti,:) = sum(XY.*nucdata)/sum(nucdata);

                    % radial averages
                    this.timeTraces.radAvgNuc{cii}(ti,:) = this.radialProfile(ti).NucAvgSeg(:,cii);
                    this.timeTraces.radAvgCyt{cii}(ti,:) = this.radialProfile(ti).CytAvgSeg(:,cii);
                    this.timeTraces.radAvgRatio{cii}(ti,:) = this.radialProfile(ti).NucCytRatio(:,cii);
                end
            end
            
            
        end

        function makeRadialAvgNoSeg(this, colimg, colnucmask, colcytmask, colmargin, ti)
            
            % makeRadialAvgNoSeg(this, colimg, colnucmask, colcytmask, colmargin, ti)
            % colcytmask can be left empty
            
            if ~exist('ti','var')
                ti = 1;
            end
            
            % masks for radial averages
            [radialMaskStack, edges] = makeRadialBinningMasks(...
                this.radiusPixel, this.radiusMicron, colmargin);
            % colType is used when multiple size colonies are processed
            % together, not here
            colType = 1; 
            
            b = this.boundingBox(ti,:);
            % crop the maskstack if the colony is clipped 
            % right now only clipped left and top
            
            xmin = max(2-b(1),1);
            xmax = min(b(2), size(colimg,2)) + xmin - 1;

            ymin = max(2-b(3),1);
            ymax = min(b(4), size(colimg,1)) + ymin - 1;

            radialMaskStack{colType} = radialMaskStack{colType}(ymin:ymax,xmin:xmax,:);
%             crop = [max(2-b(1),1), min(b(2)-b(1)+1,size(radialMaskStack{colType},2)),...
%                     max(2-b(3),1), min(b(4)-b(3)+1,size(radialMaskStack{colType},1))];
%             radialMaskStack{colType} = radialMaskStack{colType}(crop(3):crop(4),crop(1):crop(2),:);
%             
            xres = this.radiusMicron / this.radiusPixel;
            if isempty(colcytmask)
                colcytmask = imdilate(colnucmask,strel('disk',round(5/xres)))-colnucmask;
            end
            
            % do radial binning
            N = size(radialMaskStack{colType},3);
            nucradavg = zeros([N this.nChannels]);
            nucradstd = zeros([N this.nChannels]);
            cytradavg = zeros([N this.nChannels]);
            cytradstd = zeros([N this.nChannels]);
            
            % store bin edges, to be reused by segmented profiles later
            this.radialProfile(ti).BinEdges = edges{colType};

            if this.nChannels > size(colimg,3) 
                nChannels = size(colimg,3);
                if ti == 1
                    warning('img is missing channel');
                end
            else
                nChannels = this.nChannels;
            end

            for ri = 1:N
                % for some reason linear indexing is faster than binary
                colnucbinmask = find(radialMaskStack{colType}(:,:,ri) & colnucmask);
                colcytbinmask = find(radialMaskStack{colType}(:,:,ri) & colcytmask);
                
                % we only want to measure if there is something there
                npix = sum(sum((colnucmask | colcytmask) & radialMaskStack{colType}(:,:,ri)));
                npixTot = sum(sum(radialMaskStack{colType}(:,:,ri)));
                
                if npix/npixTot > 0.1

                    for ci = 1:nChannels

                        imc = colimg(:,:,ci);
                        % most primitive background subtraction: minimal value
                        % within the colony
                        % min(imc(colmaskClean)) doubles the computatation time
                        imc = imc - min(imc(:));

                        imcbin = imc(colnucbinmask);
                        nucradavg(ri,ci) = mean(imcbin);
                        nucradstd(ri,ci) = std(double(imcbin));

                        imcbin = imc(colcytbinmask);
                        cytradavg(ri,ci) = mean(imcbin);
                        cytradstd(ri,ci) = std(double(imcbin));
                    end
%                 else
%                     fprintf('x');
                end
            end
            
            this.radialProfile(ti).NucAvg = nucradavg;
            this.radialProfile(ti).NucStd = nucradstd;
            this.radialProfile(ti).CytAvg = cytradavg;
            this.radialProfile(ti).CytStd = cytradstd;
        end

        function calculateMoments(this, colimg)%, colnucmask)
            
            siz = size(colimg);
            C = round(siz(1:2)/2);
            [X,Y] = meshgrid(1:siz(1),1:siz(2));
            X = X - C(1); Y = Y - C(2);

            this.CM = {};
            this.I = {};
            
            for ci = 1:this.nChannels
                
                % for testing
                %X0 = 200; Y0 = 300; R0 = 20;
                %M = (X-X0).^2 + (X-X0).*(Y-Y0) + (Y-Y0).^2 < R0^2;

                % background intensity outside colony
                R = round(this.radiusPixel);
                mask = X.^2 + Y.^2 < R.^2;
                
                M = colimg(:,:,ci);
                bg = mean(M(~mask));
                M = double(M).*mask;

                % weight for center of mass
                Mp = M - bg;
                W = double(Mp)./sum(Mp(:));
                this.CM{ci} = [sum(sum(W.*X)) sum(sum(W.*Y))];

                % inertia tensor
                Xp = X - this.CM{ci}(1);
                Yp = Y - this.CM{ci}(2);
                this.I{ci} = [   sum(sum(W.*Xp.*Xp)) sum(sum(W.*Yp.*Xp));
                        sum(sum(W.*Yp.*Xp)) sum(sum(W.*Yp.*Yp))];
                %[V,D] = eig(I);
            end
        end
        
        % getter for dependent properties
        %---------------------------------
        
        function radius = get.radius(this)
            % for compatibility with the CellTracker.colony
            radius = this.radiusPixel;
        end
        
        % setters
        %---------------------------------
        
        function setRadius(this, radiusMicron, resolution)
            
            this.radiusMicron = radiusMicron;
            this.radiusPixel = round(radiusMicron/resolution);
        end
        
        function setID(this, ID)
            
            this.ID = ID;
        end
    end
end