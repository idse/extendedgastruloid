clear; clc; close all

scriptPath = fileparts(matlab.desktop.editor.getActiveFilename);
%[parentpath,~,~] = fileparts(scriptPath);
%[repopath,~,~] = fileparts(parentpath);
addpath(genpath(scriptPath));

dataDir = scriptPath;

%% setup
%script for stitching image grids from the Nikon or Dragonfly microscope
%stitching can be done while image acquisition is taking place, or all at
%once afterwards

bare = '_p%.4d_w%.4d_t%.4d.tif';
mipdir = fullfile(dataDir, 'MIP'); %folder in which to save MIPs
if ~exist(mipdir,'dir'), mkdir(mipdir); end

previewdir = fullfile(dataDir,'previews'); %folder in which to save previews
if ~exist(previewdir,'dir'), mkdir(previewdir); end

%set the extension; possible extensions: .nd2, .ims, .lif, .tif
ext = '.nd2';
%amount of time in seconds to wait for updates if stitching is being done
%while images are taken and transferred
nseconds = 60;
%number of colonies/stitched grids; set this manually for now, or leave
%empty to determine automatically
ngrids = [4];

gridsize = [4 6]; %size of the stitched image grid (height by width)
ppg = gridsize(1)*gridsize(2); %positions per grid
montageOverlap = 10; %percent overlap between adjacent images
montageFusionGrid = false;

%stitching method; options: weightedsample, distweight, intensityweight, randomsample
%default and generally best option is weightedsample
stitchmode = 'weightedsample';

%description of conditions in the order in which they are imaged
conditions = {'B50'};
ncond = length(conditions);
% %names of channels are determined automatically from the folder name
% [~,subdir,~] = fileparts(dataDir);
% I = regexp(subdir,'_RD[1234567890]+_');
% s = strsplit(subdir(I+1:end),'_');
% channelLabel = ['DAPI', s(2:end)];

%alternatively, channel labels could be specified manually as in:
channelLabel = {'DAPI','FOXC2','NANOG','TBX6'};

%nuclear marker channel (e.g., DAPI, HOECHST; indexed from 0)
nucChannel = 0;
disp(strcat("channels:",cell2mat(strcat(" ",channelLabel))))
disp(strcat("nuclear channel: ",channelLabel{nucChannel+1}))

mipext = '.jpg';

%% automatic settings
%if the number of grids/colonies is empty, determine it automatically
if isempty(ngrids)
    if strcmp(ext,'.nd2') || strcmp(ext,'.lif')
        listing = dir(fullfile(dataDir,['*',ext]));
        r = bfGetReader(fullfile(dataDir,listing(1).name));
        %do we also need to account for multiple files with multiple series?
        npos = r.getSeriesCount;
        seriesflag = true;
        
        if npos == 1
            npos = length(listing);
            seriesflag = false;
        end
    else
        npos = length(listing);
        seriesflag = false;
    end
    ngrids = npos/ppg;
else
    seriesflag = false;
end

%determine order of images in a stitched grid
order = NaN(ppg,2);
idx = 1;
if ~montageFusionGrid
    for ii = 1:gridsize(1)
        if mod(ii,2) == 0 %alternate direction for odd and even rows
            roworder = gridsize(2):-1:1;
        else
            roworder = 1:gridsize(2);
        end
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
else
    for ii = 1:gridsize(1)
        roworder = 1:gridsize(2);
        order(idx:idx + gridsize(2) - 1,:) = [ii*ones(gridsize(2),1),roworder(:)];
        idx = idx + gridsize(2);
    end
end

grid = NaN(gridsize);
for idx = 1:ppg
    ii = order(idx,1); jj = order(idx,2);
    grid(ii,jj) = idx;
end
disp('grid layout')
disp(grid)

fprintf("number of colonies = %d\n",ngrids)


%% do stitching
pidx = 1;
while pidx <= ngrids
    %find all the image (nd2) files in the directory
    listing = dir(fullfile(dataDir,['*',ext]));
    nfiles = numel(listing);
    fprintf('%d of %d files found\n',nfiles,ngrids*ppg)
    %if there are at least enough images to stitch the pidxth grid, stitch
    %that one and increment pidx
    if nfiles >= pidx*ppg || seriesflag
        if pidx == 1
            %assume the number of time points is the same across all files
            %(for this type of data in general it should be = 1)
            r = bfGetReader(fullfile(dataDir,listing(1).name));
            nt = r.getSizeT;
        end
        fprintf('Stitching grid #%d\n',pidx)
        for ti = 1:nt
            fprintf('t = %d\n',ti)
            imgs = cell(gridsize);
            nucmips = cell(gridsize);
            for idx = 1:ppg
                %overall index of the (idx)th image in the (pidx)th grid
                fidx = (pidx - 1)*ppg + idx;
                fprintf('.')
                
                if seriesflag
                    %if all positions are saved as series in a single file
                    r.setSeries(fidx-1);
                else
                    fname = fullfile(dataDir,listing(fidx).name);
                    r = bfGetReader(fname);
                end
                
                if idx == 1
                    m = r.getSizeY; n = r.getSizeX; nt = r.getSizeT; nz = r.getSizeZ; nc = r.getSizeC;
                    pixelOverlapY = round(m*montageOverlap/100); pixelOverlapX = round(n*montageOverlap/100);
                    stitchedSize = round([gridsize(1) - montageOverlap*(gridsize(1)-1)/100,...
                        gridsize(2) - montageOverlap*(gridsize(2)-1)/100].*[m,n]);
                end
                
                ii = order(idx,1); jj = order(idx,2);
                img = zeros(m,n,nz,nc,'uint16');
                for ci = 1:nc
                    for zi = 1:nz
                        im = bfGetPlane(r, r.getIndex(zi-1,ci-1,ti-1)+1);
                        img(:,:,zi,ci) = im;
                    end
                end
                imgs{ii,jj} = img;
                nucmips{ii,jj} = max(squeeze(img(:,:,:,nucChannel+1)),[],3);
            end
            fprintf('\n')

            if ti == 1 %double check the overlap values in pixels on the first time point
                %check pixelOverlapY
                [~,I] = max((cellfun(@(x) mean(x(end-pixelOverlapY+1:end,:),'all'),nucmips(1:end-1,:))),[],'all','linear');
                [ii,jj] = ind2sub(gridsize - [1 0],I);
                im1 = nucmips{ii,jj}(end-pixelOverlapY+1:end,:);
                im2 = nucmips{ii+1,jj}(1:pixelOverlapY,:);
                [offset, ~, ~, ~] = xcorr2fft(im1, im2);

                if abs(offset/pixelOverlapY) > 0.25
                    warning('detected pixel overlap very different from nominal value, sticking with nominal value');
                else
                    pixelOverlapY = pixelOverlapY + offset;
                end

                %check pixelOverlapX
                [~,I] = max((cellfun(@(x) mean(x(:,end-pixelOverlapX+1:end),'all'),nucmips(:,1:end-1))),[],'all','linear');
                [ii,jj] = ind2sub(gridsize - [0 1],I);
                im1 = nucmips{ii,jj}(:,end-pixelOverlapX+1:end);
                im2 = nucmips{ii,jj+1}(:,1:pixelOverlapX);
                [~, offset, ~, ~] = xcorr2fft(im1, im2);
                
                if abs(offset/pixelOverlapX) > 0.25
                    warning('detected pixel overlap very different from nominal value, sticking with nominal value');
                else
                    pixelOverlapX = pixelOverlapX + offset;
                end
            end

            upperleft = registerImageGrid_v3(nucmips, [pixelOverlapY pixelOverlapX]);
            if strcmp(stitchmode,'weightedsample')
                zidxs = randomIndicesDistWeighted(upperleft, nucmips);
            end
            
            for ci = 1:nc
                suffix = sprintf(bare,pidx-1,ci-1,ti-1);
                fname = fullfile(dataDir,['stitched',suffix]);
                mipname = fullfile(mipdir,['stitched_MIP',suffix]);
                mipname = strrep(mipname,'.tif',mipext); %save the MIP as a png
                previewname = fullfile(previewdir,['stitched_preview',suffix(1:end-4),'.jpg']);

                stitched = zeros(stitchedSize(1),stitchedSize(2),nz,'uint16');
                for zi = 1:nz
                    if zi == 1
                        mode = 'overwrite';
                    else
                        mode = 'append';
                    end
                    ims = cellfun(@(x) x(:,:,zi,ci),imgs,'UniformOutput',false);
                    
                    if strcmp(stitchmode,'distweight')
                        [stitch, ~] = stitchImageGridDistWeighted(upperleft, ims);
                    elseif strcmp(stitchmode,'intensityweight')
                        [stitch, ~] = stitchImageGridIntensityWeighted(upperleft, ims);
                    elseif strcmp(stitchmode,'randomsample')
                        [stitch, ~] = stitchImageGridRandSample(upperleft, ims);
                    elseif strcmp(stitchmode,'weightedsample')
                        stitch = applyRandomIndices(upperleft, ims, zidxs);
                    end
                    
                    %trim to projected size
                    tmp = zeros(stitchedSize,'uint16');
                    yrange = 1:min(stitchedSize(1),size(stitch,1));
                    xrange = 1:min(stitchedSize(2),size(stitch,2));
                    tmp(yrange, xrange) = stitch(yrange, xrange);
                    stitched(:,:,zi) = tmp;

                    imwrite(tmp,fname,'WriteMode',mode)
                end
                %write mip and mipidx
                [MIP,~] = max(stitched,[],3);
                if strcmp(mipext,'.jpg')
                    imwrite(im2double(MIP),mipname,'Quality',99)
                else
                    imwrite(MIP,mipname)
                end
                
                %write preview image
                imwrite(im2double(imadjust(MIP,stitchedlim(MIP))),previewname)
            end
        end
        pidx = pidx + 1;
    end
    %pause if there are not already enough files to stitch the next grid
    if numel(listing) < ngrids*ppg && numel(listing) < pidx*ppg && ~seriesflag
        disp('paused')
        timestamp=datetime('now');
        disp(timestamp)
        pause(nseconds)
    end
end

%populate and save metadata
%'nWells',           ncond,...
manualMeta = struct(...
    'posPerCondition',  ngrids/ncond,...
    'nPositions',       ngrids,...
    'montageGridSize',  gridsize,...
    'montageOverlap',   montageOverlap,...
    'xSize',            stitchedSize(2),...
    'ySize',            stitchedSize(1),...
    'nucChannel',       nucChannel,...
    'channelLabel',     {channelLabel},...
    'conditions',       {conditions});

meta = Metadata(dataDir, manualMeta);
save(fullfile(dataDir,'meta.mat'),'meta');

%% make overviews -> set the grid size and order

nc = length(channelLabel); %number of channels
%ncond = length(conditions); %number of conditions
ppc = ngrids/ncond; %positions per condition
%ppc = 20/5; %positions per condition

%set the grid size ([width height]) for displaying overviews 
gridsize = [4 3];
%order in which to display images (left to right then to the next row or
%top to bottom then to the next column)
gridorder = 'leftright';

if strcmp(gridorder,'topbottom')
    layout = [kron((1:gridsize(2))',ones(gridsize(1),1)),...
        repmat((1:gridsize(1))',gridsize(2),1)];
elseif strcmp(gridorder,'leftright')
    layout = [repmat((1:gridsize(1))',gridsize(2),1),...
        kron((1:gridsize(2))',ones(gridsize(1),1))];
end

%% load images
imgs = cell(gridsize(1),gridsize(2),nc);

for pidx = 1:ngrids
    ii = layout(pidx,1); jj = layout(pidx,2);
    if ~isnan(ii)
        for ci = 1:nc
            fprintf('.')
            fname = fullfile(mipdir,['stitched_MIP',sprintf(bare,pidx-1,ci-1,0)]);
            fname = strrep(fname,'.tif',mipext);
            im = imread(fname);
            imgs{ii,jj,ci} = im;
        end
    end
    fprintf('\n')
end
I = find(cellfun(@(x) ~isempty(x), imgs(:,:,1)),1);
[row,col] = ind2sub(size(imgs,[1,2]),I);
mn = size(imgs{row,col,1}); m = mn(1); n = mn(2);
imclass = class(imgs{row,col,1});

%% show images, separate channels
cfs = 0.075; %relative font size (font height/image height)
margin = 0.01; %margin/spacing between subplots
close all
for ci = 1:nc
    %determine contrast limits in each channel based on the brightest image
    %in that channel
    means = cellfun(@(x) mean(x,'all'),imgs(:,:,ci));
    [~,I] = max(means,[],'all','linear'); [row,col] = ind2sub(size(means),I);
    lim = stitchedlim(imgs{row,col,ci});
    figure('WindowState','maximized')
    for pidx = 1:ngrids
        cidx = floor((pidx-1)/ppc) + 1;
        ii = layout(pidx,1); jj = layout(pidx,2);
        if ~isnan(ii)
            idx = sub2ind(size(imgs,[2,1]),jj,ii);
            subplot_tight(size(imgs,1),size(imgs,2),idx,margin)
            imshow(imadjust(imgs{ii,jj,ci},lim),...
                    'InitialMagnification','fit','Border','tight')
            cleanSubplot
            %channel label
            text(50, m*(1-0.5*cfs), channelLabel{ci},...
                'Color','w','FontUnits','normalized','FontSize',cfs,...
                'FontWeight','bold','Interpreter','none')
            %condition label
            text(50, m*0.5*cfs, conditions{cidx},...
                'Color','w','FontUnits','normalized','FontSize',cfs,...
                'FontWeight','bold','Interpreter','none')
        end
    end
    drawnow
    exportgraphics(gcf, fullfile(dataDir,'previews',['overview_',channelLabel{ci},'.png']),'Resolution',300)
end
close all

%% combined colors
cfs = 0.075; %relative font size (font height/image height)
margin = 0.01; %margin/spacing between subplots
%red, green, blue, gray -> set any channel to zero if it is not used
order = [3 2 4 1];
colors = {'red','green','blue','white'};
% order = [2 1];
% colors = {'red','white'};

%combined channel label
clabel = '';
for ii = 1:length(order)
    if order(ii) > 0
        clabel = strcat(clabel,"\color{",colors{ii},"}",channelLabel{order(ii)}," ");
    end
end

%get contrast limits
lim = zeros(nc,2);
for ci = 1:nc
    means = cellfun(@(x) mean(x,'all'),imgs(:,:,ci));
    [~,I] = max(means,[],'all','linear'); [row,col] = ind2sub(size(means),I);
    lim(ci,:) = stitchedlim(imgs{row,col,ci});
end

figure('WindowState','maximized')
for pidx = 1:ngrids
    cidx = floor((pidx-1)/ppc) + 1;
    ii = layout(pidx,1); jj = layout(pidx,2);
    if ~isnan(ii)
        idx = sub2ind(size(imgs,[2,1]),jj,ii);
        subplot_tight(size(imgs,1),size(imgs,2),idx,margin)
        RGB = zeros(m,n,3,imclass);
        for ri = 1:length(order) %make RGB image
            ci = order(ri);
            if ci > 0
                im = imadjust(imgs{ii,jj,ci},lim(ci,:));
                if ri == 4
                    RGB = RGB + 0.25*repmat(im,1,1,3);
                else
                    RGB(:,:,ri) = RGB(:,:,ri) + 0.75*im;
                end
            end
        end
        if order(4) == 0
            RGB = 4/3*RGB; %scale up the brightness if there is no white channel
        end

        imshow(RGB,'InitialMagnification','fit','Border','tight')
        cleanSubplot

        %condition label
        text(50, m*0.5*cfs, conditions{cidx},...
            'Color','w','FontUnits','normalized','FontSize',cfs,...
            'FontWeight','bold','Interpreter','none')

        %channel label
        if sum(imgs{ii,jj,1},'all')>0
            ypos = m*(1-0.5*cfs); xpos = 0.01*n;
            text(xpos, ypos, clabel,...
                'FontUnits','normalized','FontSize',cfs,'FontWeight','bold')
        end
    end
end

name = 'overview_c';
for ri = 1:length(order)
    if order(ri) > 0
        name = [name,'_',channelLabel{order(ri)}]; %#ok<AGROW>
    end
end
name = [name,'.png'];
exportgraphics(gcf, fullfile(dataDir,'previews',name),'Resolution',300)
close all

