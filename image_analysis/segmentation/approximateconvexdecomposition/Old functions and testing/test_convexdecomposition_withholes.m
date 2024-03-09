clear; close all;

dataDir = 'Y:\Bohan\200319_BMPlive5by5_18h_63x\MIP';
fname = 'testmask.tif';
imname = 'stitched_MIP_p0001_w0000_t0099.tif';
fname = fullfile(dataDir,fname);
imname = fullfile(dataDir,imname);

pattern = 'stitched_MIP_p0001_w0000_t%04d.tif';
seg = imread(fname);
img = im2double(imadjust(imread(imname)));

xl = [2596 3532];
yl = [1126 2062];
seg = seg(yl(1):yl(2),xl(1):xl(2));
img = img(yl(1):yl(2),xl(1):xl(2));

% imshowpair(seg,imadjust(img))
% for ti = 1:233
%     imgname = fullfile(dataDir,sprintf(pattern,ti-1));
%     img = imread(imgname);
%     imshowpair(imadjust(img),seg)
%     title(num2str(ti))
%     pause(0.05)
% end

opts.cleanupOptions = struct('separateFused', false,...
    'clearBorder',true,...
    'minAreaStd', 1,...
    'minSolidity',0.95,...
    'minArea', 700,...
    'openSize', 5,...
    'fillholes', false);
%do initial cleanup
newseg = nuclearCleanup(seg,opts.cleanupOptions);
% imshow(double(cat(3,seg,newseg,newseg)));
imshow(cat(3,img*0.7+0.3*seg,img*0.7+0.3*newseg,img))
pause 
close all
seg = newseg;

%%
close all
tau1 = 1.3;
tau2 = 7;
s = 1;
t = 1;
opts = struct(...
    'flag',     true,... 
    'img',      img',...
    'flag3',    false);

% newseg = separate_merged_withholes(seg, tau1, tau2, s, t, opts);
newseg = separate_fused(seg, tau1, tau2, opts);

imshowpair(seg,newseg)

%%
close all
[B, L, n, A] = bwboundaries(seg);
opts = struct(...
    'flag',     true,... 
    'img',      img',...
    'flag3',    false);
[hole, parent] = find(A);
s = 1;
t = 1;

tic
[components, cuts] =...
    approx_decomp_withholes_v2(B{parent}, B(hole), s, t, opts);
toc

%% Load the large micropattern image
close all
dataDir = 'Z:\Heemskerk_Lab\Shared-03\191112_micropattern_forSeth\MIP';
readname = 'stitched_MIP_p0000_w0000_t0000_Simple Segmentation.h5';
readname = fullfile(dataDir,readname);
imname = 'stitched_MIP_p0000_w0000_t0000.tif';
imname = fullfile(dataDir, imname);
rawimg = imread(imname);
img = im2double(imadjust(imread(imname)));
seg = h5read(readname,'/exported_data');
seg = squeeze(seg) == 1;
newseg = bwareaopen(seg,3);

newseg = ~bwareaopen(~newseg, 3);
imshowpair(seg,newseg)
pause
close all
seg = newseg';

opts = struct(...
    'flag',     true,...  
    'img',      img',...
    'flag3',    false);
tau1 = 1;
tau2 = 8;

%%%%%%%
%get boundaries and hole dependencies
disp('getting boundaries')
tic
[B, ~, ~, A] = bwboundaries(seg);
toc

n = length(B);
goods = ones(n,1);

%simplify boundary contours and prevent self-intersections
disp('simplifying boundaries')
tic
for bi = 1:n
    b = B{bi};
    if size(unique(b,'rows'),1) > 3 && ~pointsAreCollinear(b)
        B{bi} = correctContours(b);
    else
        goods(bi) = 0;
    end
end
toc
goods = goods == 1;

sizes = cellfun(@(x) size(x,1), B);
[~, idx] = max(sizes);
coords = B{idx};
hidxs = A(:,idx) > 0;
hidxs = hidxs & goods;
holes = B(hidxs);

f = figure;
imshow(img')
line(coords(:,1),coords(:,2),'Color','b')
for bi = 1:length(holes)
    line(holes{bi}(:,1),holes{bi}(:,2),'Color','r')
end
pause
close(f)
%%%%%%%
%%
tic
[components, cuts, ~, ~, ~] = ...
    make_cuts_v2(coords, holes, [], [], [], tau1, tau2, opts);
toc

%% Decompose the micropattern with new function
close all; clc;
% imsize = size(seg);
% scale = 1;
% % [B, ~, ~, A] = bwboundaries(imresize(seg,imsize*scale));
% % B = cellfun(@(x) x/scale, B, 'UniformOutput', false);
% [B, ~, ~, A] = bwboundaries(seg);
% tic
% for bi = 1:length(B)
%     b = B{bi};
%     if size(unique(b,'rows'),1) > 3 && ~pointsAreCollinear(b)
%         B{bi} = correctContours(b);
%     end
% end
% toc

% sizes = cellfun(@(x) size(x,1),B);
% [~, parent] = max(sizes);
% coords = B{parent};
% children = B((A(:,parent) == 1));

tau1 = 1.2;
tau2 = 2;
s = 1;
t = 1;
opts = struct(...
    'flag',     false,...  
    'img',      img',...
    'flag3',    false);

close all
c1 = clock;
newseg = separate_fused(seg, tau1, tau2, opts);
c2 = clock;

elapsed = (c2(end-1) - c1(end-1))*60 + c2(end) - c1(end);
disp(strcat("Total elapsed time is ", num2str(elapsed), " seconds"))


imshowpair(seg,newseg)

%% visualize results
close all
CC1 = bwconncomp(seg);
idxlists = CC1.PixelIdxList;
ncomps = length(idxlists);
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
colors = distinguishable_colors(ncomps,'k');
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end

beforeimg = 0.5*cat(3,im1,im2,im3) + 0.5*cat(3,img,img,img);

CC2 = bwconncomp(newseg);
idxlists = CC2.PixelIdxList;
ncomps = length(idxlists);
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
colors = distinguishable_colors(ncomps,'k');
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end

afterimg = 0.5*cat(3,im1,im2,im3) + 0.5*cat(3,img,img,img);

subplot(1,2,1)
imshow(beforeimg)
title('Before resolution')
subplot(1,2,2)
imshow(afterimg)
title('After resolution')
set(gcf,'WindowState','maximized')

%% Try decomposing new seg again with old function
close all
tic
[newerseg, ~] = separate_merged_v3(newseg, tau1, tau2, s, t, opts);
toc

CC3 = bwconncomp(newerseg);
idxlists = CC3.PixelIdxList;
ncomps = length(idxlists);
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
colors = distinguishable_colors(ncomps,'k');
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end

aftererimg = 0.5*cat(3,im1,im2,im3) + 0.5*cat(3,img,img,img);

imshow(aftererimg)




