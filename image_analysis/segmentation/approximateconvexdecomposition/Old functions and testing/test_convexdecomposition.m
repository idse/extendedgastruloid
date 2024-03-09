clear; close all;
%select data
% dataDir = "Y:\Seth\191208_short_tracking_test\MIP";
% fname = "stitched_MIP_p0000_w0000_Simple Segmentation.h5";
% imname = "Separate times\stitched_MIP_p0000_w0000_t0001.tif";
% segidx = 1;

%other data sets to test on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dataDir = 'Y:\Seth\0_doseresponse\170726_BMPDoseResponse\MIP';
% fname = '170721_BMPDoseResponse_MIP_p0000_w0001_Simple Segmentation.h5';
% imname = '170721_BMPDoseResponse_MIP_p0000_w0001.tif';
% segidx = 2;

dataDir = 'Y:\Bohan\200319_BMPlive5by5_18h_63x\MIP';
fname = 'stitched_MIP_p0001_w0000_Simple Segmentation.h5';
imname = 'stitched_MIP_p0001_w0000_t0099.tif';
segidx = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%read image and segmentation
imname = fullfile(dataDir,imname);
fname = fullfile(dataDir,fname);
img = im2double(imadjust(imread(imname,1)));
seg = h5read(fname, '/exported_data');
seg = squeeze(seg);
seg = seg(:,:,100);
seg = seg == segidx; seg = seg';

% just for troubleshooting %%%%%%%%%
% img = img(287:1317,1754:2784);
% seg = seg(287:1317,1754:2784);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts.cleanupOptions = struct('separateFused', false,...
    'clearBorder',true,...
    'minAreaStd', 1,...
    'minSolidity',0.95,...
    'minArea', 500,...
    'openSize', 5);
%do initial cleanup
newseg = nuclearCleanup(seg,opts.cleanupOptions);
opts.cleanupOptions.separateFused = true;
seg2 = nuclearCleanup(seg,opts.cleanupOptions);
imshow(double(cat(3,seg,newseg,newseg)))
pause
close all
seg = newseg;

props = regionprops(seg,...
    'MajorAxisLength', 'MinorAxisLength','PixelIdxList');

%set parameter values
%tau1 is the threshold relative to cut length
tau1 = 1.3;
%tau2 is the absolute concavity threshold
tau2 = 15;
%increasing s biases towards making the cut to more highly concave points;
%increasing t biases towards shorter cuts
s = 1;
t = 1;
%if flag is set to true, then each cell can be examined individually; this
%is useful for determining good values of tau1 and tau2 for the given
%segmentatation
opts = struct(...
    'flag',     false,... 
    'img',      img',...
    'flag3',    false);
%do the decomposition
close all
tic
% newseg = separate_merged(seg, tau1, tau2, s, t, opts);
% [newseg2,xls,yls] = separate_merged_v2(seg, tau1, tau2, s, t, opts);
[newseg, diagnostics] = separate_merged_v3(seg, tau1, tau2, s, t, opts);
toc

%%
%visualize results side by side
CC1 = bwconncomp(seg,4);
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

CC2 = bwconncomp(newseg,4);
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

CC2 = bwconncomp(newseg2,4);
idxlists = CC2.PixelIdxList;
ncomps = length(idxlists);
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
colors = distinguishable_colors(ncomps,'k');
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end
afterimg2 = 0.5*cat(3,im1,im2,im3) + 0.5*cat(3,img,img,img);

CC2 = bwconncomp(seg2,4);
idxlists = CC2.PixelIdxList;
ncomps = length(idxlists);
im1 = zeros(size(seg)); im2 = im1; im3 = im1;
colors = distinguishable_colors(ncomps,'k');
for li = 1:ncomps
    im1(idxlists{li}) = colors(li,1);
    im2(idxlists{li}) = colors(li,2);
    im3(idxlists{li}) = colors(li,3);
end
original_separate = 0.5*cat(3,im1,im2,im3) + 0.5*cat(3,img,img,img);

xl = [1735 3242];
yl = [400 1907];

subplot(2,2,1)
imshow(beforeimg)
xlim(xl); ylim(yl);
title('Before resolution')
subplot(2,2,2)
imshow(original_separate)
xlim(xl); ylim(yl);
title('Older separate fused nuclei function')
subplot(2,2,3)
imshow(afterimg)
xlim(xl); ylim(yl);
title('Convex decomposition')
subplot(2,2,4)
imshow(afterimg2)
xlim(xl); ylim(yl);
title('Convex decomposition with curvature measure')

%% Visualize cells with different cuts
close all
saveDir = 'G:\My Drive\Research\Heemskerk lab\Data\200526_convex_decomposition';
pattern = 'example_%04d.tif';
f = figure;
list = dir(saveDir);
names = {list.name};
names = names(cellfun(@(x) contains(x,'.tif'),names));
idx = numel(names)+1;
imgs = {beforeimg,original_separate,afterimg,afterimg2};
titles = {'Before resolution','Older separate fused nuclei function',...
    'Convex decomposition','Convex decomposition with curvature measure'};
liststring = {'Yes','No'};

for xi = 1:size(xls,1)
    clf
    for si = 1:4
        subplot(2,2,si)
        imshow(imgs{si})
        ylim(xls(xi,:))
        xlim(yls(xi,:))
        title(titles{si})
    end
    set(f,'WindowState','maximized')
    [indx,~] = listdlg('ListString',{'Yes','No','Exit'},...
                'PromptString','Yes or no');
    if indx == 1
        savename = fullfile(saveDir,sprintf(pattern,idx));
        saveas(gcf,savename)
        idx = idx + 1;
    elseif indx == 3
        break
    end
end






