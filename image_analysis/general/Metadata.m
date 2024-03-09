classdef Metadata
    % a class to create a uniform metadata interface

    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties

        filenameFormat
        fileNames
        
        %tPerFile
        
        xres                % x-resolution
        yres                % y-resolution
        zres                % z-resolution
        xSize
        ySize

        nZslices

        nChannels           % number of channels
        channelNames        % cell array of channel names 
                            % when read by bioformats: detector name or dye
        excitationWavelength 

        channelLabel        % names given to channels by hand, 
                            % e.g. labeled protein

        nTime               % number of time points
        timeInterval        
        nucChannel          % nuclear label channel (e.g. DAPI or H2B)

        nPositions
        montageOverlap      % percent overlap of montage
        montageGridSize     % grid size n by m locations
        XYZ                 % position coordinates
        loop4well           % boolean if loop 4-well

        raw                 % store unprocessed metadata if necessary
        
        % to keep track of multi-well experiments
        posPerCondition     % positions per well
        conditions          % cell array of condition labels, following plate layout
    end
    properties (Dependent)
        conditionStartPos
        nWells              % number of wells/conditions
    end
   
    methods
        
        function this = Metadata(dataDir, manualMeta)
            % Metadata(dataDir, manualMeta) 
            % 
            % dataDir:      directory containing the raw data
            % manualMeta:   struct with fields matching class properties:
            % - nWells      
            % - loop4well   
            % - conditions  
            % - nucChannel  
            % - channelLabel
            %
            % and one parameter:
            % - plateLayoutNames % if files are names A1.. etc rather than a fixed format
            
            if nargin == 1 || nargin == 2
                
                % find first image file in directory
                extension = {'*.ims','*.lif','*.nd2','*tif','*.vsi'};
                
                if exist(dataDir,'file')==2 % file
                    filename = dataDir;
                    
                elseif exist(dataDir,'file')==7 % directory
                    i = 1; 
                    filelist = [];
                    while isempty(filelist) && i < numel(extension)+1
                        filelist = dir(fullfile(dataDir,extension{i}));
                        i = i+1;
                    end
                    if ~isempty(filelist)
                        filename = fullfile(dataDir,filelist(1).name);
                    else
                        error('no image files found in data directory');
                    end
                end

                % read metadata from image file (and possibly corresponding
                % txt file)
                this = this.read(filename);
                
                % set format of raw data image files
                [C, matches] = strsplit(filename,'_[Ffwptm][0-9]+','DelimiterType','RegularExpression','CollapseDelimiters',false);
                if ~isempty(matches)
                    
                    nmatches = numel(matches);
                    this.filenameFormat = C{1};
                    for mi = 1:nmatches
                        ndigits = numel(matches{mi}) - 2;
                        matchformat = [matches{mi}(1:2) '%.' num2str(ndigits) 'd'];
                        this.filenameFormat = [this.filenameFormat matchformat C{mi+1}];
                    end
                else
                    this.filenameFormat = filename;
                end
            end
            % if there is manual metadata:
            if nargin == 2
                
                % copy over manual metadata
                fn = fieldnames(manualMeta);
                for i = 1:numel(fn)
                   if isprop(this, fn{i})
                       this.(fn{i}) = manualMeta.(fn{i});
                   end
                end
%                 if ~isfield(manualMeta,'nWells')
%                     warning('please provide nWells, assuming nWells=1');
%                     manualMeta.nWells = 1;
%                 end
                if ~isfield(manualMeta,'plateLayoutNames')
                    manualMeta.plateLayoutNames = false;
                end
                if ~isfield(manualMeta,'loop4well')
                    this.loop4well = false;
                end

                % if there are individual filenames rather than a format,
                % we use a list this.fileNames and clear the filenameFormat
                if manualMeta.plateLayoutNames
                    fileNames = getPlateFileNames(filelist, this.conditions);
                    this.fileNames = fileNames;
                    this.nPositions = sum(cellfun(@numel, this.fileNames(:)));
                    if isempty(this.posPerCondition)
                        this.posPerCondition = this.nPositions / this.nWells;
                    end
                    this.filenameFormat = []; % to prevent confusion
                else
                    % could generate a cell array of filenames from
                    % filenameFormat here, seems unnecessary right now
                    if isempty(this.posPerCondition)
                        this.posPerCondition = this.nPositions / this.nWells;
                    end
                    if this.nPositions ~= this.posPerCondition*this.nWells
                        warning('positions ~= posPerConditon*nWells');
                    end
                    this.fileNames = cell(size(this.conditions));
                    for wi = 1:this.nWells
                        
                        for cpi = 1:this.posPerCondition
                            
                            pi = this.posPerCondition*(wi-1) + cpi - 1; % count from zero
                            this.fileNames{wi}{cpi} = sprintf(this.filenameFormat, pi);
                        end
                    end
                end
            end
        end
        
        function this = read(this, filename)
            % read metadata from file using bioformats
            %
            % read(filename)
            
            [dataDir,barefname,extension] = fileparts(filename);
            
            r = bfGetReader(filename);
            omeMeta = r.getMetadataStore();
            
            this.xSize = omeMeta.getPixelsSizeX(0).getValue();
            this.ySize = omeMeta.getPixelsSizeY(0).getValue();
            
            this.nZslices = r.getSizeZ();
            this.nTime = r.getSizeT();
            
            this.xres = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM));
            this.yres = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM));
            if this.nZslices > 1
                this.zres = double(omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM));
            end
            
            dt = omeMeta.getPixelsTimeIncrement(0);
            if ~isempty(dt)
                this.timeInterval = double(dt.value);
            else
                this.timeInterval = [];
            end
            
            this.nChannels = r.getSizeC();
            this.channelNames = {};
            for ci = 1:this.nChannels 
                this.channelNames{ci} = char(omeMeta.getChannelName(0,ci-1));
            end
            
            this.excitationWavelength = {};
            for ci = 1:this.nChannels 
                lambda = omeMeta.getChannelExcitationWavelength(0,ci-1);
                if ~isempty(lambda)
                    this.excitationWavelength{ci} = round(10^3*double(lambda.value(ome.units.UNITS.MICROM)));
                end
            end

            %omeMeta.getPixelsType(0);
            %temporary fix, don't think we actually use meta.raw for much
            try
                this.raw = char(omeMeta.dumpXML());
            catch
                this.raw = [];
            end
            
            
            % ---- SPECIAL FOR LASX -------------
            nSeries = r.getSeriesCount();
            if nSeries > 1 && strcmp(extension,'.lif')
                this.nPositions = nSeries;
                warning('multiple series in file, assuming nPositions = nSeries');
            end
            
            % ---- SPECIAL FOR FUSION ---------------
            fusionMetaFilename = fullfile(dataDir,[barefname '_metadata.txt']);
            
            if exist(fusionMetaFilename,'file')
            
                disp('reading Fusion metadata file');
                
                fid = fopen(fusionMetaFilename);
                if fid == -1
                    error('file not found');
                end
                tline = '';
                i = 1;

                while ~isnumeric(tline) 
                    if contains(tline,'[FieldDefinition]') 
                        xline = strsplit(fgets(fid),'X=');
                        yline = strsplit(fgets(fid),'Y=');
                        zline = strsplit(fgets(fid),'ReferenceZ=');
                        this.XYZ(i,:) = [str2num(xline{2}) str2num(yline{2}) str2num(zline{2})];
                        i = i+1;
                    end
                    tline = fgets(fid);
                end
                this.XYZ = this.XYZ*10^3; %mm->um
                this.nPositions = size(this.XYZ, 1);
                fclose(fid);
            else
                warning('Fusion meta file not found');
            end
        end
        
        function save(this)
            % save this object to a mat file
            %
            % save()
            %
            % e.g. raw filename is 1.oib -> stores metadata in same place
            % under 1_metadata.mat

            [datadir,barefname] = fileparts(this.filename);
            metafname = fullfile(datadir,[barefname '_metadata']);
            save(metafname, 'this');
        end
        
        function displayPositions(this)
            % display positions 
            %
            % displayPositions();
            %
            % positions are center of field of view?
            
            XYZ = this.XYZ;

            %scatter(XYZ(:,1), XYZ(:,2))

            for i = 1:size(XYZ,1)
                text(XYZ(i,1), XYZ(i,2),num2str(i))
                w = 1024*this.xres;
                h = 1024*this.yres;
                rectangle('Position',[XYZ(i,1)-w/2,XYZ(i,2)-h/2,w,h])
            end
            axis([min(XYZ(:,1))-w max(XYZ(:,1))+w min(XYZ(:,2))-h max(XYZ(:,2))+h])
            axis equal
            axis off

            % XYZmean = mean(XYZ);
            % hold on
            % scatter(XYZmean(:,1), XYZmean(:,2),'r')
            % hold off
        end
        
        function conditionStartPos = get.conditionStartPos(this)
            % positions at which conditions start
            conditionStartPos = 1 + [0 cumsum(this.posPerCondition)];
            conditionStartPos = conditionStartPos(1:end-1);
        end
        
        function nWells = get.nWells(this)
            nWells = numel(this.conditions);
        end
    end
end