% VERSION HISTORY
%
% 231115 : Idse - change to deal with self-stitched files
%

clear; close all;

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);

% raw data location, modify if not the same as the location of this script
dataDir = scriptPath;
cd(dataDir);

% define channels and conditions
%--------------------------------------
manMeta = struct();
manMeta.nucChannel = 0;
manMeta.channelLabel = {'DAPI','ISL1','pSMAD1','FOXF1'};
manMeta.conditions = {'BMP','X','LDN'};
manMeta.nChannels = numel(manMeta.channelLabel);

% proper preprocessing should make it unnecessary to define the resolution
manMeta.zres = 2;
manMeta.xres = 0.353;
manMeta.yres = 0.353;
manMeta.nTime = 1;
%manMeta.nZslices = 17;

% define micropatterned 'MP' or 'Disordered'
imageType = 'MP'; 
% radii of micropatterns in micron for each condition (ignored if not 'MP')
radii = [700 700 700]/2;  
% margin in microns to keep around the nominal radius (for masking cells)
Rmarg=55; 

% different cases for input files
%--------------------------------------

% self-stitched - channels are split
if ~isempty(dir(fullfile(dataDir,'stitched*.tif')))
    filelistAll = {};
    filenameFormat = "stitched_p%.4d_w%.4d_t%.4d.tif";
    xsectfname_prefix_format = 'stitched_p%.4d';
    
% Dragonfly stitched - channels are not split
elseif ~isempty(dir(fullfile(dataDir,'*FusionStitcher*.ims')))
    filelistAll = dir(fullfile(dataDir,'*FusionStitcher*.ims'));
    % FINISH
end

% define files per condition
%--------------------------------------

filelist = {};

% case 1: fixed number of positions per condtion 
posPerCondition = 4;
manMeta.nWells = numel(manMeta.conditions);
manMeta.posPerCondition = repmat(posPerCondition, [1 manMeta.nWells]);
manMeta.nPositions = posPerCondition*manMeta.nWells;
conditionStartPos = 1 + [0 cumsum(manMeta.posPerCondition)];
conditionStartPos = conditionStartPos(1:end-1);

if ~isempty(filelistAll)
    for i = 1:numel(manMeta.conditions)
        filelist{i} = filelistAll((i-1)*posPerCondition+1:i*posPerCondition);
    end
% else
%     filelist = {};
% 
%     for condi = manMeta.nWells
%        filelist{condi} = []; 
%     end
% 
%     for pi = 1:manMeta.nPositions
%         condi = find(pi >= conditionStartPos,1,'last');
%         filelist{condi} = [filelist{condi} sprintf(filenameFormat, pi)];
%     end
end

% % case 2: define files per condition manually 
% % required if variable number of positions per conditions
% filelist{1} = dir(fullfile(dataDir,'*80K*.ims'));
% filelist{2} = dir(fullfile(dataDir,'*160K*.ims'));
% filelist{3} = dir(fullfile(dataDir,'*240K*.ims'));
% filelistAll = cat(2,filelist(:));
% 
% % define number of positions for each condition
% manMeta.posPerCondition = cellfun(@numel,filelist);

% starting position of each conditioni is meta.conditionStartPos

% reload previous metadata or create new
%--------------------------------------
metafile = fullfile(dataDir,'meta.mat');
if exist(metafile)
    load(metafile,'meta');
    disp(['loading meta file ' metafile]);
else
    meta = Metadata(dataDir, manMeta);
    save(metafile, 'meta');
    disp(['created meta file ' metafile]);
end

% reload previous results
%--------------------------------------
posfile = fullfile(dataDir,'positions.mat');
if exist(posfile)
    load(posfile, 'positions');
    disp(['loading positions file ' posfile]);
end

statsfile = fullfile(dataDir,'stats.mat');
if exist(statsfile)
    load(statsfile, 'stats');
    disp(['loading stats file ' statsfile]);
end

countfile = fullfile(dataDir,'counts.mat');
if exist(countfile)
    load(countfile, 'counts');
    disp(['loading count file ' countfile]);
end


% ratio between z and xy resolution
scalefactor = round(meta.zres / meta.xres);


%% cleanup stack for specific positions

zcorrection = true;
zcorrfactor = [1.02 1 1 1];

tic
%parpool(3)
%posidx = 5;%[1,5,9];
posidx = [1,8];
parfor i = 1:numel(posidx)
%for i = 1:numel(posidx)
    fname = char(sprintf(filenameFormat, posidx(i)));
    cleanupStack(dataDir, fname, zcorrection, zcorrfactor)
end
toc

%micropatternPieVisConditions(dataDir, positions, options, meta);


%% make xsections of different channels

%pList = [1,2,4,9];
pList = [9];
RGBset = [4 3];

% parameters for cross-section
%tol = [0.01 0.97; 0.01 0.995; 0.38 0.92; 0.01 0.98];
tol = [0.01 0.97; 0.01 0.995; 0.01 0.999; 0.01 0.999];
%yinds = {1350:1360, 1320:1330, 1270:1280, 1230:1240, 1200:1210, 1400:1410, 1430:1440,1490:1500}; %change yinds if preferred
%yinds = {1370:1380, 1390:1400, 1290:1300, 1260:1270, 1400:1410, 1450:1460,1510:1520};
yinds = {1350:1360};
% pList = 1;
% yinds = {1350:1360};
% RGBset = [2 3 4];

count = 0; % use first position to define lim

for pi = pList
    count = count + 1;

    % self-stitched (channels split)
    if filenameFormat.startsWith('stitched')

        i = 0;
        for ci = 0:meta.nChannels-1
            i = i + 1;

            fullfname = fullfile(dataDir,sprintf(filenameFormat, pi-1, ci, 0));
            [fpath,barefname,ext] = fileparts(fullfname);
            barefname = char(barefname);
            filteredfname = fullfile(fpath, 'filtered', ['stitched_filtered' barefname(9:end) char(ext)]);

            if exist(filteredfname)
                disp('using filtered stack');
                fullfname = filteredfname;
            else
                warning('no filtered stack');
            end
            
            stack{i} = readStack(fullfname);
            % readStack mistakes z for t in stitched tiffs
            stack{i} = permute(stack{i}, [1 2 3 5 4]); 
        end
        img = cat(3,stack{:});
        xsectfname_prefix = sprintf(xsectfname_prefix_format, pi-1);
        
    % not self-stitched 
    else
        fname = filelistAll{pi}.name;
        img = readStack(fullfile(dataDir,fname));
        xsectfname_prefix = fname(1:end-4);
    end

    for di = 1:2
        if di == 1
            direction = 'x';
        else
            direction = 'y';
        end
        for i = 1:numel(yinds)

            options = struct('tol',tol, 'index', yinds{i}, 'RGBset', RGBset,...
                'color', 'GM','direction',direction);

            if count == 1
                disp('setting Ilim');
                options.Ilim = makeXsection(img, xsectfname_prefix, meta, options);
            else
                makeXsection(img, xsectfname_prefix, meta, options);
            end
            close;
        end
    end
end

%%
%--------------------------------------------------------------------------
% INITIAL QC
%--------------------------------------------------------------------------

% make xsections of segmentation 

pList = 1:4;
yinds = {1350:1360};

tic
for pi = pList %1:manMeta.nWells
    
    ci = meta.nucChannel;
    
    % load image data
    fullfname = fullfile(dataDir,sprintf(filenameFormat, pi-1, ci, 0));
    [fpath,barefname,ext] = fileparts(fullfname);
    barefname = char(barefname);
    filteredfname = fullfile(fpath, 'filtered', ['stitched_filtered' barefname(9:end) char(ext)]);

    if exist(filteredfname)
        disp('using filtered stack');
        fullfname = filteredfname;
    else
        warning('no filtered stack');
    end

    img = readStack(fullfname);

    % load segmentation
    fullfname = char(fullfile(dataDir,sprintf(filenameFormat, pi-1, ci, 0)));
    prefix = fullfname(1:end-4);
    masksfname = [prefix '_masks.mat'];
    masks = load(masksfname);
    label = uint16(masks.bgmask*0);
    for i = 1:numel(masks.masks)
        label(masks.masks(i).nucmask)=i;
    end
    
    %change yinds and xinds if preferred
    i = 1;
    xinds = 1:size(img,2);
    ymarg=0;
    yind = yinds{i}(1)-ymarg:yinds{i}(end)+ymarg;
    res = makeSegmentationCrossSection(img, label, meta, yind, xinds);
    %saveas(gcf, [prefix  'segmentationXsection_y' num2str(yind(1)) ' ' num2str(yind(end)) '.png']);
    close;
    imwrite(res.segoverlay, [sprintf(xsectfname_prefix_format, pi-1)  '_segmentation_' num2str(yind(1)) ' ' num2str(yind(end)) '_yz.png']);
    figure()
    MIP = max(squeeze(img),[],3);
    imshow(imadjust(MIP,stitchedlim(MIP)));
    y = [min(yinds{i}),max(yinds{i})];
    x = [min(xinds),max(xinds)];
    line(x,y,'Color','r','LineWidth',2)
    saveas(gcf, [sprintf(xsectfname_prefix_format, pi-1)  'segmentationXsection_overview_' num2str(yind(1)) ' ' num2str(yind(end)) '_yz.png']);
    close;
end
toc


%% load previously saved XZ xsection and scatter cell centers to check segmentation

figure,

pi = 9; % position 
i = 1; % cross-section index

% ci = meta.nucChannel; % channel
%
fname = [sprintf(xsectfname_prefix_format, pi-1) '_' meta.channelLabel{ci+1} '_' num2str(yinds{i}([1 end])) '_yz.png'];
fname = [sprintf(xsectfname_prefix_format, pi-1)  '_segmentation_' num2str(yind(1)) ' ' num2str(yind(end)) '_yz.png'];
xsect = imread(fullfile(dataDir,fname));

% in y range 
ymarg = 13; % a margin is needed to get all the positions (why?)
inyrange = (positions(pi).cellData.XY(:,2) < yinds{i}(end) + ymarg) & (positions(pi).cellData.XY(:,2) > yinds{i}(1) - ymarg);
sum(inyrange)

% z is flipped in the xsection image so we have to flip the coordinate too
Z = size(xsect,1) - positions(pi).cellData.Z(inyrange)*scalefactor; 
X = positions(pi).cellData.XY(inyrange,1);

s = 100;
imshow(xsect)
%imshow(1-imadjust(mat2gray(xsect)))

hold on 
scatter(X,Z,s,'r.');
hold off

%% overlay quantification for QC on XZ cross-section

pi = 9;
ci = 2;
condi = find(pi >= conditionStartPos,1,'last');
% ci = find(pi >= conditionStartPos,1,'last'); ???
fname = [sprintf(xsectfname_prefix_format, pi-1) '_' meta.channelLabel{ci+1} '_' num2str(yinds{i}([1 end])) '_yz.png'];
xsect = imread(fullfile(dataDir,fname));

inyrange = (positions(pi).cellData.XY(:,2) < yinds{i}(end) + ymarg) & (positions(pi).cellData.XY(:,2) > yinds{i}(1) - ymarg);
Z = size(xsect,1) - positions(pi).cellData.Z(inyrange)*scalefactor; 
X = positions(pi).cellData.XY(inyrange,1);

C = positions(pi).cellData.nucLevel(:,ci+1);

tol = 0.01;
n = length(C);
Cs = sort(C);
Cmin = Cs(max(ceil(n*tol),1));
Cmax = Cs(floor(n*(1-tol)));
C(C < Cmin) = Cmin; 
C(C > Cmax) = Cmax;
C = mat2gray(C);

imshow(mat2gray(xsect))
s = 300;
hold on 
scatter(X,Z,s,C(inyrange),'.');
hold off
%saveas(gcf, fullfile(dataDir,['quantificationOverlay_' meta.channelLabel{ci} '_' meta.conditions{condi} '.png']));

%% overlay quantification for QC on XY MIP

piSerious = 1;%[1,5,9];
ci = 4; % define nucChannel you want to check
checkOPT = 'nuc'; % choose from 'nuc' 'cyto' 'NCratio'

for pi = piSerious
    condi = find(pi >= conditionStartPos,1,'last');

    tol = 0.01;
    img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);

    % saturate ALL PIXEL according to tol
    imgGREY = imadjust(img,stretchlim(img,[tol,1-tol]),[]);
    imgBLACK = imbinarize(img, 1);

    nucLevel = positions(pi).cellData.nucLevel;
    cytLevel = positions(pi).cellData.cytLevel;
    background = positions(pi).cellData.background;
    XY = positions(pi).cellData.XY;

    cNuc = nucLevel(:,ci)-background(ci);
    cCyt = cytLevel(:,ci)-background(ci);
    NC = cNuc./cCyt;

    if strcmp(checkOPT,'nuc')
        z = cNuc;
    elseif strcmp(checkOPT,'cyto')
        z = cCyt;
    elseif strcmp(checkOPT,'NCratio')
        z = NC;
    end

    % saturate ALL CELLS according to tol
    n = length(z);
    zs = sort(z);
    zmin = zs(max(ceil(n*tol),1));
    zmax = zs(floor(n*(1-tol)));
    z(z < zmin) = zmin; z(z > zmax) = zmax;
    z = (z-min(z))/(max(z)-min(z));

    figure()
    imshow(imgGREY)
    hold on
    s = scatter(XY(:,1), XY(:,2), 50, z*[1, 1, 1], "filled");
    hold off
    %s.MarkerEdgeColor = 'm';
    title(['Intensity in ',meta.channelLabel{ci}])
    saveas(gcf,['Intensity_check_' meta.channelLabel{ci} '_' manMeta.conditions{condi} '_Pos' num2str(pi) '.jpg']);
end

%% check if scatter plots of markers in different samples in same conditions look similar

condi = 1; % define the condition you want to check
i = 2; % define two channels you want to compare
j = 3; % define two channels you want to compare

X = stats.nucLevel{condi}(:,i);
Y = stats.nucLevel{condi}(:,j);
X = log(1+X/mean(X));
Y = log(1+Y/mean(Y));
scatter(X, Y, 10, stats.sample{condi}, 'filled')
% xlim([0 1.5]);
% ylim([0 3]);

%%
%--------------------------------------------------------------------------
% QUANTIFY MARKERS AND DENSITY IN 3D (MP only)
%--------------------------------------------------------------------------

%% load previously saved xsection and scatter cell centers in polar coordinates

pi = 1;
ci = 0;
i = 1;

fname = [sprintf(xsectfname_prefix_format, pi-1) '_' meta.channelLabel{ci+1} '_' num2str(yinds{i}([1 end])) '_yz.png'];
xsect = imread(fullfile(dataDir,fname));
imshow(xsect)

% in y range 
scalefactor = round(meta.zres / meta.xres);

% z is flipped in the xsection image so we have to flip the coordinate too
Z = size(xsect,1) - positions(pi).cellData.Z*scalefactor; 
R = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));

hold on 
s = 10;
scatter(positions(pi).center(1) + R,Z,s,'r.');
scatter(positions(pi).center(1) - R,Z,s,'r.');
hold off

%% scatter colored for colony per condition

clf 
hold on
colors = lines(6);
s = 10;
condi = 3;

i = 1;
for pi = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1

    Z = size(xsect,1) - positions(pi).cellData.Z*scalefactor; 
    R = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));

    scatter(+ R,size(xsect,1)-Z,s,colors(i,:),'.');
    scatter(- R,size(xsect,1)-Z,s,colors(i,:),'.');
    i = i+1;
    disp([num2str(pi) ' ncells: ' num2str(numel(R))])
end
legend({'1','2','3','4'})
    
hold off
axis equal
ylim([0 240]);
xlim([-max(R+Rmarg) max(R+Rmarg)]);
cleanSubplot;
saveas(gcf, fullfile(dataDir,['radialScatter_colonies_' meta.conditions{condi} '.png']));

%% scatter colored for marker

condi = 2;

for ci = 2:4
    clf 
    hold on

    for pi = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1

        Z = size(xsect,1) - positions(pi).cellData.Z*scalefactor; 
        R = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));
        C = positions(pi).cellData.nucLevel(:,ci);

        tol = 0.01;
        n = length(C);
        Cs = sort(C);
        Cmin = Cs(max(ceil(n*tol),1));
        Cmax = Cs(floor(n*(1-tol)));
        C(C < Cmin) = Cmin; 
        C(C > Cmax) = Cmax;

        scatter(+ R,size(xsect,1)-Z,s,C,'.');
        scatter(- R,size(xsect,1)-Z,s,C,'.');
    end
    
    hold off
    axis equal
    ylim([0 200]);
    xlim([-max(R+Rmarg) max(R+Rmarg)]);
    cleanSubplot;
    saveas(gcf, fullfile(dataDir,['radialScatter_' meta.channelLabel{ci} '_' meta.conditions{condi} '.png']));
end

%% scatter cells in YR colored for density

%figure,
clf

pi = 10;
Z = positions(pi).cellData.Z*scalefactor; 
R = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));

[f,xi] = ksdensity([R,Z],'Support','positive','Bandwidth',0.15);
np = sqrt(size(xi,1));
Xi = reshape(xi(:,1),[np,np]); 
Zi = reshape(xi(:,2),[np,np]); 
F = reshape(f,[np,np]);
density = interp2(Xi,Zi,F,R,Z);
%density = density./R;

tol = 0.01;
n = length(density);
densitys = sort(density);
densitymin = densitys(max(ceil(n*tol),1));
densitymax = densitys(floor(n*(1-tol)));
density(density < densitymin) = densitymin; 
density(density > densitymax) = densitymax;

% hold on
% scatter(positions(pi).center(1) + R,Z,15,density,'filled');
% scatter(positions(pi).center(1) - R,Z,15,density,'filled');
% hold off
% cleanSubplot
% axis equal
% ylim([min(Z), max(Z)])
% xlim([positions(pi).center(1) - max(R), positions(pi).center(1) + max(R)])

hold on
scatter(R,Z,15,density,'filled');
hold off
cleanSubplot
axis equal
ylim([min(Z), max(Z)])
xlim([0, max(R+Rmarg)])
yticks([0 max(Z)]);


%% plot density of cells in cross section

clf
%figure,
pi = 10;
Z = positions(pi).cellData.Z*scalefactor; 
R = sqrt(sum((positions(pi).cellData.XY - positions(pi).center).^2,2));

[f,xi] = ksdensity([R,Z],'Support','positive','Bandwidth',0.15);
np = sqrt(size(xi,1));
Ri = reshape(xi(:,1),[np,np]); 
Zi = reshape(xi(:,2),[np,np]); 
F = reshape(f,[np,np]);

[Rg,Zg] = meshgrid(1:round(max(R)),1:round(max(Z)));

density = interp2(Ri,Zi,F,Rg,Zg);
%density = density./Rg;
imshow(flipud(density),[]);
colormap turbo;

hold on 
s = 20;
scatter(R,round(max(Z))-Z, s, '.k');
hold off

%% plot c

C = positions(pi).cellData.nucLevel(:,3);
F = scatteredInterpolant(R,Z,C);
expression = F(Rg,Zg);

imshow(expression,[])



%% CLEAN UP positions & VISUALIZE cells

checkPos = true; % set if you want to VISUALIZE cells
posCLEANUP = false; % set if you want to do clean up

badidxPos = cell(numel(positions),1);
badidxCon = cell(numel(manMeta.conditions),1);
for pos= 1:numel(positions)
    
    condi = ceil(pos/meta.posPerCondition);
    i = 1; % define which channel you do CLEAN UP according to

    background = positions(pi).cellData.background;
    nucLevel = positions(pos).cellData.nucLevel(:,i)-background(i);
    nucArea = positions(pos).cellData.nucArea(:,1);

    stdFOLD = 2;
    % throwout brightest DAPI (dying/dead cells)
    badidxPos{pos} = (nucLevel > mean(nucLevel) + stdFOLD*std(nucLevel));
    % throw out darkest DAPI (usually segmentation error)
    badidxPos{pos} = badidxPos{pos} | (nucLevel < mean(nucLevel) - stdFOLD*std(nucLevel));
    % throw out cells with too small nucArea, also can use nucAxis, nucCircularity etc
    badidxPos{pos} = badidxPos{pos} | (nucArea < mean(nucArea) - stdFOLD*std(nucArea));
    if strcmp(imageType,'MP')
        % throw out cells outside the colony edges
        badidxPos{pos} = badidxPos{pos}; (sqrt(sum((positions(pos).cellData.XY-positions(pos).center).^2,2))*meta.xres > positions(pos).radiusMicron);
    end
    if posCLEANUP == 0
        badidxCon{condi} = [badidxCon{condi};badidxPos{pos}];
    end

    if checkPos == 1 % check how CLEAN UP works
        figure()
        img = max(positions(pos).loadImage(dataDir,i-1,1),[],3);
        img = mat2gray(img, [150 4000]);
        imshow(img)
        hold on
        XY = positions(pos).cellData.XY;
        XYbad = positions(pos).cellData.XY(badidxPos{pos},:);
        scatter(XY(:,1),XY(:,2), 7,'filled','c')
        scatter(XYbad(:,1),XYbad(:,2), 7,'filled','m')
        hold off
        saveas(gcf, fullfile(dataDir, ['CLEAN UP position' num2str(pos) '.png']));
    end

    if posCLEANUP == 1 % throw out badIdx cells
        positions(pos).ncells = positions(pos).ncells - sum(badidxPos{pos});
        positions(pos).cellData.cytMaskArea = positions(pos).cellData.cytMaskArea(~badidxPos{pos},:);
        positions(pos).cellData.cytLevel = positions(pos).cellData.cytLevel(~badidxPos{pos},:);
        positions(pos).cellData.NCratio = positions(pos).cellData.NCratio(~badidxPos{pos},:);
        positions(pos).cellData.XY = positions(pos).cellData.XY(~badidxPos{pos},:);
        positions(pos).cellData.nucArea = positions(pos).cellData.nucArea(~badidxPos{pos},:);
        positions(pos).cellData.nucOrientation = positions(pos).cellData.nucOrientation(~badidxPos{pos},:);
        positions(pos).cellData.nucMajorAxis = positions(pos).cellData.nucMajorAxis(~badidxPos{pos},:);
        positions(pos).cellData.nucMinorAxis = positions(pos).cellData.nucMinorAxis(~badidxPos{pos},:);
        positions(pos).cellData.nucCircularity = positions(pos).cellData.nucCircularity(~badidxPos{pos},:);
        positions(pos).cellData.nucZ = positions(pos).cellData.nucZ(~badidxPos{pos},:);
        positions(pos).cellData.nucLevel = positions(pos).cellData.nucLevel(~badidxPos{pos},:);
        badidxCon{condi} = [badidxCon{condi};badidxPos{pos}(~badidxPos{pos},:)];
    end

end

%% get stats

stats = cellStats(positions, meta, positions(1).dataChannels);
% automaticly get thresholds
confidence = 0.95;
conditions = 1:numel(meta.conditions);
whichthreshold = []; %[1 1 2]
stats.getThresholds(confidence, conditions, whichthreshold);

% Make nucHistogram channel by channel
conditionIdx = 1:numel(meta.conditions); % [1,3,5] Specify conditions for Histograms
tolerance = 0.01;
nbins = 50;
stats.makeHistograms(nbins, tolerance);
for channelIdx = 1:numel(manMeta.channelLabel)
    options = struct('channelIndex',channelIdx, 'cytoplasmic', false,...
        'cumulative',false, 'time', 1,...
        'conditionIdx',conditionIdx,...
        'titlestr',meta.channelLabel{channelIdx}...
        );
    figure()
    stats.plotDistributionComparison(options)
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{channelIdx(1)} '_dist.png']));
end
%ylim([0 0.02])
%xlim([0 1000])

save(statsfile, 'stats');
%stats.exportCSV(meta); export to CSV

%% Create combined nucHistogram channel by channel, If mannually thresholding, use this as a referrence

nucLevelCombined = [];
XYCombined = [];
sampleIdx = [];
for i = 1:numel(manMeta.conditions)
    nucLevel = stats.nucLevel{i};
    nucLevelCombined = [nucLevelCombined; nucLevel];
end
for i=1:numel(manMeta.channelLabel)
    figure()
    histogram(nucLevelCombined(:,i));
    xlabel('Intensity');
    ylabel('Population');
    title([manMeta.channelLabel{i},' All'])
    %     xlim([0,3000]);
    %     ylim([0,3000]);
    saveas(gcf, fullfile(dataDir, [meta.channelLabel{i} '_tot_dist.jpg']));
end

%% make intensity radial profile (MP Only)

% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', method 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
options = struct('nucChannels', [1 2 3], ...'ratioMode', 'N:C',...
    'normalize', true, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
    'log1p',false,'edgeDistance',false ...
    );
% options = struct('edgeDistance',false);
% options = struct('ratioChannels', 1, 'ratioMode', 'N:C',...
%     'normalize', true, 'std', true, 'FontSize', 15, 'legend',true,...
%     'pointsPerBin', 200, 'interpmethod', 'nearest', 'stdMethod', 'perColony',...
%     'log1p',false);

res = {};
% first round without normalization (to collect overall max and min of all), keep normalize = false

conditionindices = [1,2,3];%[13 16:18];
%conditionindices = [13:15];

for condi = conditionindices%1:numel(meta.conditions)
    clf
    
    tempoptions = options;
    tempoptions.normalize = false;
    tempoptions.conditionIdx = condi;
    tempoptions.colonyRadius = radii(condi);
    
    res{condi} = plotRadialProfiles(stats, meta, tempoptions);
 
    min(res{condi}.ratio_profile)

    if condi == conditionindices(1)
        maxvalsNucAll = max(res{condi}.nuc_profile);
        minvalsNucAll = min(res{condi}.nuc_profile);

        maxvalsCytAll = max(res{condi}.cyt_profile);
        minvalsCytAll = min(res{condi}.cyt_profile);

        maxvalsRatioAll = max(res{condi}.ratio_profile);
        minvalsRatioAll = min(res{condi}.ratio_profile);
    else
        maxvalsNucAll = max(max(res{condi}.nuc_profile), maxvalsNucAll);
        minvalsNucAll = min(min(res{condi}.nuc_profile), minvalsNucAll);

        maxvalsCytAll = max(max(res{condi}.cyt_profile), maxvalsCytAll);
        minvalsCytAll = min(min(res{condi}.cyt_profile), minvalsCytAll);

        maxvalsRatioAll = max(max(res{condi}.ratio_profile), maxvalsRatioAll);
        minvalsRatioAll = min(min(res{condi}.ratio_profile), minvalsRatioAll);
    end
end
close all

% second round with normalization
nuclimits = [minvalsNucAll' maxvalsNucAll'];
options.nuclimits = nuclimits;

cytlimits = [minvalsCytAll' maxvalsCytAll'];
options.cytlimits = cytlimits;

ratiolimits = [minvalsRatioAll' maxvalsRatioAll'];
options.ratiolimits = ratiolimits;

for condi = conditionindices
    options.conditionIdx = condi;
    options.colonyRadius = radii(condi);
    figure()
    res{condi} = plotRadialProfiles(stats, meta, options);
    ylim([0 1.1]);
    %ylim([-0.0195 -0.016]);
    xlim([0 radii(condi)]);
    axis square
    title(meta.conditions{condi})
    saveas(gcf, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.png'])
    writetable(res{condi}.datatable, ['radialProfile_' meta.channelLabel{:} '_' meta.conditions{condi} '.csv'])
    %close;
end

%% radial intensity - conditions combined

markerChannels = options.nucChannels+1;
for channel = [3]%markerChannels
    figure()
    %colors = lines(numel(conditionindices));
    colors = [0.8500    0.3250    0.0980;; 0.4660    0.6740    0.1880; 0    0.4470    0.7410; 0 0 0]
    lw = 4;
    fs = 25;
    legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.nuc_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        %legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
        legendstr = [legendstr,{[meta.conditions{condi}]}];
    end
    i = 0;
    for condi = conditionindices
        i = i+1;
        
        plerr = res{condi}.nuc_profile(:,channel) + res{condi}.nuc_profile_colstd(:,channel);
        minerr = res{condi}.nuc_profile(:,channel) - res{condi}.nuc_profile_colstd(:,channel);
        % plerr = res{condi}.nuc_profile(:,channel) + res{condi}.nuc_profile_std(:,channel);
        % minerr = res{condi}.nuc_profile(:,channel) - res{condi}.nuc_profile_std(:,channel);
        good = ~isnan(plerr);
        fill([res{condi}.r(good),fliplr(res{condi}.r(good))],[plerr(good)', fliplr(minerr(good)')],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    hold off
    xlim([0 radii(condi)]);
    xlim([0 400]);
    ylim([0 1.25]);
    xlabel('r (\mum)')
    %xlabel('edge distance (\mum)')
    if options.log1p
        ylabel('log( intensity )');
    else
        ylabel('intensity');
    end
    legend(p,legendstr,'FontSize',20,'Location','northwest')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    axis square;
    box off;
    cleanSubplot(fs,2);
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.png']));
    export_fig(fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.pdf']));
end

%% radial intensity - conditions combined for ratio

markerChannels = options.ratioChannels+1;

for channel = markerChannels
    figure()
    colors = lines(numel(conditionindices));

    lw = 4;
    fs = 25;
    legendstr = [];
    p = [];
    hold on
    i = 0;
    for condi = conditionindices
        i = i+1;

        plot(res{condi}.r, res{condi}.ratio_profile(:,channel),'LineWidth',lw,'Color',colors(i,:));
        legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
    end
    i = 0;
    for condi = conditionindices
        i = i+1;
        
        plerr = res{condi}.ratio_profile(:,channel) + res{condi}.ratio_profile_std(:,channel);
        minerr = res{condi}.ratio_profile(:,channel) - res{condi}.ratio_profile_std(:,channel);
        good = ~isnan(plerr);
        fill([res{condi}.r(good),fliplr(res{condi}.r(good))],[plerr(good)', fliplr(minerr(good)')],...
            colors(i,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    hold off
    xlim([0 400]);
    ylim([-0.1 1.5]);
    xlabel('r (\mum)')
    if options.log1p
        ylabel('log( intensity )');
    else
        ylabel('N:C intensity');
    end
    legend(p,legendstr,'FontSize',15,'Location','northwest')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    axis square;
    saveas(gcf, fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.png']));
    export_fig(fullfile(dataDir, ['radialintensity_combine_', meta.channelLabel{channel} '.pdf']));
end

%% make intensity radial profile, one colony (MP Only)

% 1. set 'nucChannels' 'cytChannels' 'ratioChannels' to plot different value
% 2. default interpmethod is 'linear', 'nearest' is for sparse colony center
% 3. stdMethod: 'perColony' 'neighborCells', usually perColony give less STD
% 4. decrease pointsPerBin when you analysis single colony
posList = [1:4];
norMethodPos = 1; % 1 - normalize individually; 2 - normalize according to all position
options = struct('nucChannels', 1,... %'ratioMode', 'C:N',...
    'normalize', false, 'std', true, 'FontSize', 15, 'legend',true,...
    'pointsPerBin', 100, 'interpmethod', 'nearest', 'stdMethod', 'perColony','log1p',false);
% first round without normalization (to collect overall max and min of all), keep normalize = false
count = 0;
for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition(1));
    resOP = {};
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    %metaOP.nWells = 1;
    metaOP.posPerCondition = 1;

    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);

    clf
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    close all

    count = count+1;
    if count == 1
        maxvalsNuc = max(resOP.nuc_profile);
        minvalsNuc = min(resOP.nuc_profile);

        maxvalsCyt = max(resOP.cyt_profile);
        minvalsCyt = min(resOP.cyt_profile);

        maxvalsRatio = max(resOP.ratio_profile);
        minvalsRatio = min(resOP.ratio_profile);
    else
        maxvalsNuc = max(maxvalsNuc, max(resOP.nuc_profile));
        minvalsNuc = min(minvalsNuc, min(resOP.nuc_profile));

        maxvalsCyt = max(maxvalsCyt, max(resOP.cyt_profile));
        minvalsCyt = min(minvalsCyt, min(resOP.cyt_profile));

        maxvalsRatio = max(maxvalsRatio, max(resOP.ratio_profile));
        minvalsRatio = min(minvalsRatio, min(resOP.ratio_profile));
    end
end

if norMethodPos == 2
    nuclimits = [minvalsNuc' maxvalsNuc'];
    options.nuclimits = nuclimits;
    cytlimits = [minvalsCyt' maxvalsCyt'];
    options.cytlimits = cytlimits;
    ratiolimits = [minvalsRatio' maxvalsRatio'];
    options.ratiolimits = ratiolimits;
end

options.normalize = true; % set normalization here

for posidx = posList
    condi = ceil(posidx/manMeta.posPerCondition(1));
    options.conditionIdx = 1;
    options.colonyRadius = radii(condi);
    metaOP = meta;
    metaOP.conditions = meta.conditions(condi);
    %metaOP.nWells = 1;
    metaOP.posPerCondition = 1;
    statsOP = cellStats(positions(posidx), metaOP, positions(1).dataChannels);
    figure()
    resOP = plotRadialProfiles(statsOP, metaOP, options);
    ylim([0 1.1]);
    xlim([0 radii(condi)]);
    axis square
    title([meta.conditions{condi} ' Pos' num2str(posidx)])
    saveas(gcf, ['radialProfileOP_' meta.channelLabel{:} '_' meta.conditions{condi} '_Pos' num2str(posidx) '.png'])
end


%% scatter plot subpopulation on XY MIP to check thresholds

condi = 1; % condition index
cpi = 0; % position within condition index
s = 30;

pi = cpi + meta.conditionStartPos(condi);
%subDir = filelist{condi}(pi).folder;

%channelThresholds = [0, stats.thresholds(2), stats.thresholds(3), stats.thresholds(4)]; % use automatic thresholds
channelThresholds = [0 2300 1800 1330]; % manually adjust thresholds

nucLevel = positions(pi).cellData.nucLevel;
background = positions(pi).cellData.background;
positive = {};
for i = 1:meta.nChannels
    positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
end

for ci = 3%:numel(manMeta.channelLabel)
    
    % s = strsplit(positions(pi).filename,{'_FusionStitcher','.ims'});
    % prefix = [s{:}];
    %
    % zdir = [prefix '_zslices'];
    % img = imread(fullfile(subDir, zdir, sprintf([prefix '_MIP_w%.4d.tif'], ci-1)));
    
    img = max(positions(pi).loadImage(dataDir,ci-1,1),[],3);
    img = imadjust(img,stitchedlim(img));
    %img = mat2gray(img, [150 4000]);
    [X,Y] = meshgrid(1:size(img,2),1:size(img,1));
    if strcmp(imageType,'MP')
        R = sqrt((X - positions(pi).center(1)).^2 + (Y - positions(pi).center(2)).^2);
        mask = R > positions(pi).radiusPixel + Rmarg/meta.xres;
        img(mask) = max(img(:));
    end

    figure,
    imshow(img)
    hold on
    %scatter(positions(pi).cellData.XY(:,1),positions(pi).cellData.XY(:,2),'filled')
    XY = positions(pi).cellData.XY(positive{ci},:);
    scatter(XY(:,1),XY(:,2),s,'filled','r')
    %XY = positions(pi).cellData.XY(AP2Cp,:);
    %scatter(XY(:,1),XY(:,2),'x','g')
    XY = positions(pi).cellData.XY(positive{2},:);
    %scatter(XY(:,1),XY(:,2),'x','b')
    if strcmp(imageType,'MP')
        scatter(positions(pi).center(1), positions(pi).center(2),500,'.','g')
    end
    hold off
    title(meta.channelLabel{ci})
    saveas(gcf,['subpopulation_' meta.channelLabel{ci} '_Pos' num2str(pi) '.jpg']);
end
% save new manually adjusted thresholds
for pi = 1:numel(positions)
    positions(pi).cellData.channelThresholds = channelThresholds;
end
save(fullfile(dataDir,'positions'), 'positions');
stats.thresholds = channelThresholds;
save(statsfile, 'stats');


%% scatter plot subpopulation on XZ xsection to check thresholds

i = 1;
pi = 1;
ci = 2;

fname = [sprintf(xsectfname_prefix_format, pi-1) '_' meta.channelLabel{ci+1} '_' num2str(yinds{i}([1 end])) '_xz.png'];
xsect = imread(fullfile(dataDir,fname));

inyrange = (positions(pi).cellData.XY(:,2) < yinds{i}(end) + ymarg) & (positions(pi).cellData.XY(:,2) > yinds{i}(1) - ymarg);
Z = size(xsect,1) - positions(pi).cellData.Z(inyrange)*scalefactor; 
X = positions(pi).cellData.XY(inyrange,1);

nucLevel = positions(pi).cellData.nucLevel;
background = positions(pi).cellData.background;
positive = {};
for i = 1:meta.nChannels
    positive{i} = nucLevel(:,i) - background(i) > channelThresholds(i);
end
idx = positive{ci+1}(inyrange);

imshow(mat2gray(xsect))
s = 300;
hold on 
scatter(X(idx),Z(idx),s,'.r');
hold off
%saveas(gcf, fullfile(dataDir,['quantificationOverlay_' meta.channelLabel{ci} '_' meta.conditions{condi} '.png']));

%% probability of being positive in RZ by condition (MP Only)
%----------------------------------------------------------------------

% DO WE NEED SOME SORT OF Z-ALIGNMENT?
% NOW DO THIS FOR INTENSITIES?

condi = 3;
conditionPositions = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1;
pos = positions(conditionPositions);
%options = struct('half',true);
%options = struct('half',true,'color','RGB', 'channels',[4 2 3]);

options = struct('half',true,'color','MG', 'channels',[3 4 NaN],'flip',false,'FontSize',34.5);
radpos2Dcomb = plotRadialPositive2D(pos, meta, stats.thresholds, options);
saveas(gcf, fullfile(dataDir, ['radialpositive2D_' meta.conditions{condi} '_' options.color '.png']));
%writetable(radpos2Dcomb.dataTable,fullfile(dataDir, ['radialpositive2D_' meta.conditions{condi} '.csv']))
%export_fig(fullfile(dataDir, ['radialpositive2D_' meta.conditions{condi} '_' options.color '.pdf']));%,'-nocrop');


%% overlay of mean profile for different conditions

radpos2Dcomb = {};

for condi = 1:3

    conditionPositions = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1;
    pos = positions(conditionPositions);
    options = struct('half',true);
    radpos2Dcomb{condi} = radialPositive2D(pos, meta, stats.thresholds, options);
end

figure,
%clf
fs = 14;
lw = 2;
colors = gray(4);
colors = colors(end-1:-1:1,:);
hold on
for condi = 1:3
    B = radpos2Dcomb{condi}.boundaries{1};
    if ~options.half
        B(:,2) = B(:,2) - Rmax;
    else
        B = B(B(:,1) > 1 & B(:,2) > 1,:);
    end
     plot([B(:,2)', -flipud(B(:,2))', B(1,2)]*meta.xres,...
          [B(:,1)', flipud(B(:,1))', B(1,1)]*meta.xres,'-','LineWidth',lw,'Color',colors(condi,:))
%     plot(B(:,2)'*meta.xres, B(:,1)'*meta.xres,'-','LineWidth',lw,'Color',colors(condi,:))
%     plot(-B(:,2)'*meta.xres, B(:,1)'*meta.xres,'-','LineWidth',lw,'Color',colors(condi,:))
end
hold off
cleanSubplot(fs,lw)
axis equal
axis off
ylim([0 1.1]*radpos2Dcomb{3}.Zmax*meta.xres);
xlim([-1 1]*radpos2Dcomb{3}.Rmax*meta.xres);
legend(meta.conditions,'Location','North','Orientation','horizontal');
legend('boxoff')
%saveas(gcf, fullfile(dataDir, 'shapeVsCondition.png'));


%% 2D RZ intensity profile
%----------------------------------------------------------------------

condi = 3;

conditionPositions = meta.conditionStartPos(condi):meta.conditionStartPos(condi)+meta.posPerCondition(condi)-1;
pos = positions(conditionPositions);
options = struct('half',true);
radpos2Dcomb = radialPositive2D(pos, meta, stats.thresholds, options);
%radprof2Dcomb = radialProfile2D(pos, meta, options);

radprof2Dcomb = radialProfile2D(stats, meta, condi, options);

% imshow(radprof2Dcomb.nuc_profile(:,:,1),[])
% set(gca,'YDir','normal')
% 
% hold on
% % scatter plot some points
% Npts = 10000;
% R = radprof2Dcomb.R;
% Z = radprof2Dcomb.Z;
% idx = 1+round(rand([1 Npts])*(numel(R)-1));
% %idx = 1:numel(R);
% scatter(R(idx), Z(idx),1,'.r')
% hold off

% % BOUNDARY SEEMS OFF - but actually no worse than for density estimate from
% % radialPositive2D
% imshow(radpos2Dcomb.density_scaled,[])
% set(gca,'YDir','normal')
% 
% hold on
% % scatter plot some points
% Npts = 10000;
% R = radprof2Dcomb.R;
% Z = radprof2Dcomb.Z;
% idx = 1+round(rand([1 Npts])*(numel(R)-1));
% %idx = 1:numel(R);
% scatter(R(idx), Z(idx),1,'.r')
% hold off


DAPI = mat2gray(radprof2Dcomb.nuc_profile(:,:,1));
R = mat2gray(radprof2Dcomb.nuc_profile(:,:,2));
G = mat2gray(radprof2Dcomb.nuc_profile(:,:,3));
B = mat2gray(radprof2Dcomb.nuc_profile(:,:,4));

%figure,imshow(cat(3,R,G,B))

%figure, 
RGB = 1 - cat(3, DAPI + R, DAPI + G, DAPI + B)/2;
%RGB = 1 - cat(3, DAPI + G+B, DAPI + R+B, DAPI + R+G)/2;

h = size(RGB,1);
w = size(RGB,2);
    
figure('Position',[0 0 w, h]),
    
imshow(RGB)
set(gca,'YDir','normal')

hold on
B = radpos2Dcomb.boundaries{1};
if ~options.half
    B(:,2) = B(:,2) - Rmax;
else
    B = B(B(:,1) > 1 & B(:,2) > 1,:);
end
plot(B(:,2), B(:,1),'-','LineWidth',3,'Color','k')
hold off

saveas(gcf, fullfile(dataDir, ['radialintensity2D_' meta.conditions{condi} '.png']));


%% radial probability of being positive by condition (MP Only)

Pall = {};
Pstdall = {};
xall = {};

options = struct('edgeDistance',false);

for condi = [1,2,3]%1:numel(meta.conditions)
    combos = {}; % [2 4],[2 3],[3 4]};
    stats.markerChannels = 2:4;
    [P,Pc,x,Pstd,Pcstd] = radialPositive(stats, condi, meta, combos, options);
    Pall{condi} = P;
    Pstdall{condi} = Pstd;
    xall{condi} = x;
    ylim([0 1]);
    %saveas(gcf, fullfile(dataDir, ['radialpositive_' meta.conditions{condi} '_', meta.channelLabel{2:end} '.png']));
end

%% radial probability of being positive channel by channel (MP Only)

markerChannels = [4];
for channel = markerChannels
    figure()
    %colors = lines(numel(manMeta.conditions));
    %colors = lines(8);
    %colors = colors([2 5 1 3 4 6 7 8],:);
    %colors = colors([1 2 4 3 5 6 7 8],:);
    colors = [0.8500    0.3250    0.0980; 0.4660    0.6740    0.1880; 0    0.4470    0.7410; 0 0 0 ]
    lw = 4;
    fs = 25;
    legendstr = [];
    p = [];
    for condi = [1,2,3]%1:numel(meta.conditions)
        Pspe = Pall{condi}(channel,:);
        Pstd_spe = Pstdall{condi}(channel,:);
        xspe = xall{condi};

        ptemp = plot(xspe,Pspe,'LineWidth',lw,'Color',colors(condi,:));
        %legendstr = [legendstr,{[meta.conditions{condi},' ',meta.channelLabel{channel}]}];
        legendstr = [legendstr,{[meta.conditions{condi}]}];
        p = [p,ptemp];
        hold on
        good = ~isnan(Pspe);
        fill([xspe(good),fliplr(xspe(good))],...
            [Pspe(good) + Pstd_spe(good), fliplr(Pspe(good) - Pstd_spe(good))],...
            colors(condi,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    xlim([0 max(xall{condi})]);
    ylim([0 1]);
    xlabel('r (\mum)')
    ylabel('positive fraction ');
    legend(p,legendstr,'FontSize',20,'Location','northwest')
    set(gca,'FontSize', fs)
    set(gca,'FontWeight', 'bold')
    set(gca, 'LineWidth', 2);
    axis square;
    box off;
    cleanSubplot(fs,2);
    saveas(gcf, fullfile(dataDir, ['radialpositive_combine_', meta.channelLabel{channel} '.png']));
    export_fig(fullfile(dataDir, ['radialpositive_combine_', meta.channelLabel{channel} '.pdf']));
end

%% count populations

combo = [4 3 2];
conditionsidx = 1:numel(manMeta.conditions);
counts = countPopulations(positions, meta, stats, dataDir, combo, conditionsidx);
save(countfile, 'counts');

%% make pretty scatter plot

for condi = [1]%1:numel(meta.conditions)%[1,2,3]%

    close all;
    options = struct();

    options.selectDisType = [1, 1, 1, 1]; % 1 - nuc, 2 - cyto, 3 - NCratio, 4 - CNratio
    options.channelThresholds = stats.thresholds;%[0 2300 1800 1330];%stats.thresholds; % re-define channelThres when plot cyto or NCratio
    %options.checkPos = badidxCon; % whether you want to check defined cells in scatter plot
    options.showThreshold = [true true true true];
    options.minimal = false;
    options.edgeDistance = true;
    options.imageType = imageType;
    options.conditionIdx = condi;
    options.channelCombos = {[4 3 2]};%, [4 2 3]};%, [4 2 3]}; % {[3 2 4], [4 3 2], [4 2 3]}

    options.axisLabel = meta.channelLabel;
    options.channelMax = exp([3 3 1.5 3])-1; %xlim or ylim
    options.log1p = [true true false true];
    %options.log1p = [false false false false];
    if strcmp(imageType,'MP')
        options.radiusMicron = max(xall{condi}); %positions(meta.conditionStartPos(condi)).radiusMicron;%
    end
    options.conditionsCombined = false;
    %options.showThreshold = [false false true false];

    results = scatterMicropattern(stats, meta, dataDir, options);
    %writetable(results.dataTable, fullfile(dataDir, ['scatterData_' meta.conditions{condi} '.csv']))
end


%% pretty channel combo pie (MP Only)

options = struct();
options.channels = [2 4];
options.pieOrder = [1 2];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.01 0.995; 0.01 0.995; 0.01 0.995];
options.tol = [0.1 0.995; 0.1 0.995; 0.1 0.995; 0.1 0.995];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.scalebar = false;
options.color = 'CMY';

posidx = 1:12; % the first position sets the lookup table

for pi = posidx

    condi = find(pi >= meta.conditionStartPos,1,'last');
    cpi = pi - meta.conditionStartPos(condi) + 1;
    %fname = filelist{condi}(cpi).name;
    %subDir = filelist{condi}(cpi).folder;
    %subDir = fullfile(dataDir, 'MIP');
    
    if pi == posidx(1)
        Ilim = micropatternPieVis(dataDir, positions(pi), options);
    else
        options.Ilim = Ilim;
        micropatternPieVis(dataDir, positions(pi), options);
    end
end

%% pretty channel combo pie (MP Only)

% options = struct();
% options.channels = [2 3];
% options.pieOrder = [1 2];
% % tolerances in original order of channels
% options.tol = [0.1 0.995; 0.01 0.995; 0.01 0.995; 0.01 0.995];
% %options.tol = [0.1 0.995; 0.1 0.995; 0.1 0.995; 0.1 0.995];
% options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
% options.radialMargin = 50; % margin in disk shaped mask of colony in micron
% options.scalebar = false;
% options.color = 'CYM';
% 


%pi = [2,4,9];
options = struct();
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 110; % margin in disk shaped mask of colony in micron
options.color = 'GM';
options.channels = [4 3];
options.pieOrder = [2 1];
options.scalebar = false;
posidx = [2,4,9]; % the first position sets the lookup table

%micropatternPieVis(dataDir, positions(pi), options);

for pi = posidx

    % condi = find(pi >= meta.conditionStartPos,1,'last');
    % cpi = pi - meta.conditionStartPos(condi) + 1;
    %fname = filelist{condi}(cpi).name;
    %subDir = filelist{condi}(cpi).folder;
    %subDir = fullfile(dataDir, 'MIP');

    if pi == posidx(1)
        Ilim = micropatternPieVis(dataDir, positions(pi), options);
    else
        options.Ilim = Ilim;
        micropatternPieVis(dataDir, positions(pi), options);
    end
end

% for pi = posidx
% 
%     condi = find(pi >= meta.conditionStartPos,1,'last');
%     cpi = pi - meta.conditionStartPos(condi) + 1;
%     %fname = filelist{condi}(cpi).name;
%     %subDir = filelist{condi}(cpi).folder;
%     %subDir = fullfile(dataDir, 'MIP');
% 
%     if pi == posidx(1)
%         Ilim = micropatternPieVis(dataDir, positions(pi), options);
%     else
%         options.Ilim = Ilim;
%         micropatternPieVis(dataDir, positions(pi), options);
%     end
% end


%%
pi = 1;
options = struct();
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 50; % margin in disk shaped mask of colony in micron
options.color = 'MG';
options.channels = [2 4];
options.pieOrder = [1 2];
micropatternPieVis(dataDir, positions(pi), options);


%% multi condition pie (MP Only)

options = struct();
options.channels = [4 3];
% tolerances in original order of channels
options.tol = [0.01 0.995; 0.1 0.995; 0.01 0.995; 0.01 0.995];
options.outsideMargin = 50; % margin outside the disk shaped mask in pixels
options.radialMargin = 110; % margin in disk shaped mask of colony in micron
options.color = 'GM';
options.scalebar = true;

% options.positionIdx = [1 5 9];
% options.positionIdx = [3 6 11];
options.positionIdx = [2 9];

micropatternPieVisConditions(dataDir, positions, options, meta);


%%
% ------------------------------------------------------------------------
% cell density and confluence
% ------------------------------------------------------------------------

%% density plot condition by condition (MP Only)

rAll = cell(numel(meta.conditions),1);
densityAll = cell(numel(meta.conditions),1);

numberAll = zeros(manMeta.posPerCondition(1), numel(manMeta.conditions));
for condi = 1:numel(meta.conditions)
    rAll{condi} = res{condi}.r';
    density = [];
    for pos = 1:manMeta.posPerCondition
        posidx = (condi-1)*manMeta.posPerCondition(1)+pos;
        densityTem = res{condi}.celldensity{posidx}';
        density = [density,densityTem];
        numberAll(pos,condi) = positions(posidx).ncells;
    end
    densityAll{condi} = density;
end

PI = 3.1415926;
densityColonyAll = numberAll/(PI*350^2)*10^6;
numberAllMean = mean(numberAll,1);
numberAllStd = std(numberAll,0,1);

c = lines(5);
fs = 20;
p = [];
figure()
for condi = 1:numel(meta.conditions)
    r = rAll{condi};
    densityMEAN = mean(densityAll{condi}*10^6,2); % cells/mm^2
    densitySTD = std(densityAll{condi}*10^6,0,2); % cells/mm^2
    hold on
    errorbar(r, densityMEAN, densitySTD,'--','LineWidth',0.5, 'Color', c(condi,:));
    ptemp = plot(r, densityMEAN,'-','LineWidth',3, 'Color', c(condi,:), 'DisplayName', manMeta.conditions{condi});
    hold off
    p = [p,ptemp];
end
xlim([0 350])
ylim([0 30000])
legend(p,'Location','best');
xlabel('edge distance ( um )')
ylabel('cell density (cells/mm^2)');
set(gca,'FontSize', fs)
set(gca,'FontWeight', 'bold')
set(gca, 'LineWidth', 2);
saveas(gcf, ['radialProfile_density' '.png'])

%% Calculate density & confluency of not confluent position

posNoConfluent = [1,5,9];
Area = zeros(numel(positions),1);
Confluency = ones(numel(positions),1);
Density = zeros(numel(positions),1);
for pi = posNoConfluent
    condi = ceil(pi / meta.posPerCondition);
    cpi = pi - (condi-1)*meta.posPerCondition;
    img = max(positions(pi).loadImage(dataDir, 0, 1),[],3);
    img = imadjust(img,stitchedlim(img));
    imgBW = imbinarize(img);
    imgBWCLEAN = bwareaopen(imgBW,1000); % remove small particles
    diskSize = 30; % play with disk size
    imgBW2 = imfill(imclose(imgBWCLEAN,strel('disk',diskSize)),'holes');
    figure()
    imshowpair(imgBWCLEAN,imgBW2,'montage')
    saveas(gcf,['Area_' manMeta.conditions{condi} '_Position_' num2str(pi) '.jpg']);
    Area(pi) = bwarea(imgBW2)*meta.xres*meta.yres/10^6;
end
for pi = 1:numel(positions)
    if Area(pi) == 0
        Area(pi) = stats.area*100;
    else
        Confluency(pi) = Area(pi)/(stats.area*100);
    end
    Density(pi) = positions(pi).ncells/Area(pi); % cells/mm^2
end
DensityAll = zeros(numel(manMeta.conditions),1);
DensitySTD = zeros(numel(manMeta.conditions),1);
coloniesPerSize = cellfun(@numel,filelist);
posbins = [0 cumsum(coloniesPerSize)];
for condi = 1:numel(manMeta.conditions)
    DensityAll(condi) = mean(Density(posbins(condi)+1:posbins(condi+1)));
    DensitySTD(condi) = std(Density(posbins(condi)+1:posbins(condi+1)));
end
