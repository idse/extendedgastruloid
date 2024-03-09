clear; clc; close all;

dataDir = 'Y:\Seth\Segmentation_benchmarking\images';
maskDir = 'Y:\Seth\Segmentation_benchmarking\masks';
fijiDir = 'Y:\Seth\Segmentation_benchmarking\imagej_watersheds';
ellipseDir = 'Y:\Seth\Segmentation_benchmarking\\Ellipse fitting';
list = dir(fullfile(dataDir,'*.tif'));
fnames = {list.name};

excluded = {...
    'IXMtest_A12_s7_w1EAEEA614-51ED-43B3-A4FF-088730911E4C';...
    'IXMtest_E05_s2_w15CCE97F8-3F02-412E-8BF2-FC92972DDA1A';...
    'IXMtest_F13_s7_w13C1B1D8C-293E-454F-B0FD-6C2C3F9F5173';...
    'IXMtest_F22_s6_w1F4C7ADE4-B68D-4D30-A063-722B87AA2DA1';...
    'IXMtest_I04_s9_w16A5CC270-8B92-42EE-AA4A-855776F7D46B';...
    'IXMtest_L01_s2_w1E5038251-DBA3-44D0-BC37-E43E2FC8C174';...
    'IXMtest_L10_s6_w12D12D64C-2639-4CA8-9BB4-99F92C9B7068'};

%option for whether to visualize results image-by-image
showims = false;
%number of images to be evaluated:
nimgs = length(fnames) - length(excluded);

%make structs to store benchmarking info for each method:
fields = {'FP', 'FN', 'TP', 'DSC', 'JI', 'Precision', 'Recall',...
    'Accuracy', 'Fmeasure', 'centroid', 'trueCentroid'};

methods = {'ilastik', 'distwatershed', 'ellipsefit', 'oldmethod',...
    'convexdecomp', 'combined'};

%initialize structure arrays within a larger struct
allstats = struct;
allstats.filenames = cell(nimgs,1);
for mi = 1:length(methods)
    allstats.(methods{mi})(nimgs) = struct;
end

%ilastik hysteresis separation options
opts = struct(...
    'core',     0.95,...
    'final',    0.5,...
    'minSize',  10,...
    'maxSize',  10000000);

%relative concavity threshold
tau1 = 0.7;
%tau2 is the absolute concavity threshold
tau2 = 1.5;
decompopts = struct(...
    'flag',         false,...
    'flag3',        false,...
    'diagnostics',  false);

m = 2; n = 4;

figure
set(gcf,'WindowState','maximized')
tic
for ti = 1:length(fnames)
    if ~showims
        disp(strcat("Image ", num2str(ti)))
    end
    
    %determine filenames
    fname = fnames{ti};
    name = fname(1:end-4);
    if ~any(cellfun(@(x) strcmp(name,x), excluded))
        %parse filenames
        fname = fullfile(dataDir,fname);
        segname = fullfile(dataDir, [name,'_Probabilities.h5']);
        maskname = fullfile(maskDir, [name, '.png']);
        
        %load image, ilastik segmentation, provided mask
        img = 0.5*imadjust(im2double(imread(fname)));
        seg = h5read(segname,'/exported_data');
        probs = squeeze(seg(1,:,:))';
        seg = squeeze(seg(1,:,:) > 0.5)';
        mask = imread(maskname); mask = squeeze(mask(:,:,1));
        imsize = size(seg);
        
        %make a 3D matrix to hold data for all segs (in order found in
        %methods list above)
        allsegs = zeros(imsize(1), imsize(2), length(methods));
        
        %Fiji watershed
        im3 = imread(fullfile(fijiDir,[name,'_Simple_Segmentation-dist-watershed.tif']));
        allsegs(:,:,2) = im3 > 0;

        %Ilastik hysteresis threshold
        nucseg = ilastik_watershed(probs, opts);
        allsegs(:,:,1) = nucseg;
        
        %Ellipse-fitting method (use ellipse centroids as seeds for
        %marker-controlled watershed)
        M = readmatrix(fullfile(ellipseDir, [name, '.csv']));
        xcol = 3; ycol = 4;
        x = M(:, xcol); x(x < 1) = 1; x(x > imsize(2)) = imsize(2);
        y = M(:, ycol); y(y < 1) = 1; y(y > imsize(1)) = imsize(1);
        inds = sub2ind(imsize, round(y), round(x));
        seeds = zeros(size(seg));
        seeds(inds) = 1;
        Ds = mat2gray(bwdist(seeds));
        D = bwdist(~seg);
        Dc = imcomplement(mat2gray(D));
        Dtot = Ds + Dc;
        V = imimposemin(Dtot, seeds);
        L = watershed(V);
        L(~seg) = 0;    
        eseg = L > 0;
        allsegs(:,:,3) = eseg;

        %Old erosion + watershed-based method
        nucmask = oldSeparateFusedNuclei(seg);
        allsegs(:,:,4) = nucmask;
        
        %Convex decomposition
        newseg = separate_fused(seg,tau1,tau2,decompopts);
        allsegs(:,:,5) = newseg;

        %Ilastik hysteresis + convex decomposition
        newerseg = separate_fused(nucseg, tau1, tau2, decompopts);
        allsegs(:,:,6) = newerseg;
        
        %calculate and store all benchmarking data
        allstats.filenames{ti} = name;
        for mi = 1:length(methods)
            [stats, ~, ~] = detection_accuracy(allsegs(:,:,mi), mask, struct);
            for fi = 1:length(fields)
                allstats.(methods{mi})(ti).(fields{fi}) = stats.(fields{fi});
            end
        end
        
        %optionally show results for each iteration
        if showims
            %Raw pixel-wise mask from Ilastik
            subplot(m,n,1)
            [stats, ~, ~] = detection_accuracy(seg, mask, struct);
            [im1, ~] = visualize_nuclei(seg,img);
            imshow(im1)
            title('Raw Ilastik mask')
            xlabel(strcat("Accuracy = ", num2str(stats.Accuracy)))
            
            %Ground truth segmentation
            subplot(m,n,2)
            [im2, ~] = visualize_nuclei(mask,img);
            imshow(im2)
            title('Ground truth')
            
            %visualize performance accuracy
            for mi = 1:length(methods)
                subplot(m,n,mi+2)
                cla
                imshow(allsegs(:,:,mi))
                hold on
                stats = allstats.(methods{mi})(ti);
                xyTP = stats.trueCentroid(stats.TP,:);
                scatter(xyTP(:,1),xyTP(:,2),15,'g','filled')
                xyFN = stats.trueCentroid(stats.FN,:);
                scatter(xyFN(:,1),xyFN(:,2),15,'b','filled')
                xyFP = stats.centroid(stats.FP,:);
                scatter(xyFP(:,1),xyFP(:,2),15,'r','filled')
                title(methods{mi})
                xlabel(strcat("Accuracy = ", num2str(stats.Accuracy)))
            end
            
            imtitle = name;
            imtitle(strfind(name,"_")) = " ";
            sgtitle(imtitle)
            pause
        end
        
    else
        disp(strcat("Skipped ", name))
    end
end
toc
close all

%% Visualize aggregate results
close all
groups = {'Precision', 'Recall', 'Accuracy', 'Fmeasure'};
method = {'ilastik', 'distwatershed', 'ellipsefit', 'oldmethod',...
    'convexdecomp', 'combined'};

Y = zeros(4,6);
errlow = zeros(4,6);
for ii = 1:4
    for jj = 1:6
        vals = 100*[allstats.(method{jj}).(groups{ii})];
        Y(ii,jj) = mean(vals,'omitnan');
        %Use SEM error bars:
        errlow(ii,jj) = std(vals,'omitnan')/sqrt(numel(vals));
    end
end

X = categorical(groups);
X = reordercats(X, groups);

% bar(X,Y);
b = bar(Y);
set(gca, 'XTickLabel', groups);
ylabel('Measure percentage')
title('Detection performance across methods for 193 images')

hold on
% Find the number of groups and the number of bars in each group
ngroups = size(Y, 1);
nbars = size(Y, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, Y(:,i), errlow(:,i), 'k', 'linestyle', 'none');
end

legend(b, method, 'Location', 'southeast')
hold off

%% Visualize distribution of F-measures
close all
figure
set(gcf,'WindowState','maximized')
for jj = 1:6
    subplot(2,3,jj)
    vals = 100*[allstats.(method{jj}).(groups{ii})];
    histogram(vals,50)
    ylim([0, 60])
    xlim([40, 100])
    title(method{jj})
end

sgtitle('Distribution of F-measure across images for each approach')

%% Visualize distribution of DSC
% close all
figure
set(gcf,'WindowState','maximized')
for jj = 1:6
    subplot(2,3,jj)
    vals = cell2mat(transpose({allstats.(method{jj}).DSC}));
    disp(numel(vals))
    histogram(vals)
    ylim([0 2500])
    title(method{jj})
end

sgtitle('Distribution of DSC across all nuclei for each approach')






