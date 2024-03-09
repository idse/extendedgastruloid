function newNuclearMask = oldSeparateFusedNuclei(nuclearMask, options)

if ~exist('options','var')
    options = struct();
end
if ~isfield(options,'minSolidity')
    minSolidity = 0;
else
    minSolidity = options.minSolidity;
end
if ~isfield(options,'minAreaStd')
    minAreaStd = 1;
else
    minAreaStd = options.minAreaStd;
end

CC = bwconncomp(nuclearMask);
if minSolidity > 0
    stats = regionprops(CC, 'ConvexArea', 'Area');
    convexArea = [stats.ConvexArea];
else
    stats = regionprops(CC, 'Area');
end
area = [stats.Area];

% alternative to minAreaStd that is more robus to large numbers of
% fusions but requires more thought
if isfield(options,'minFusedArea')
    minFusedArea = options.minFusedArea;
else
    minFusedArea = mean(area) + minAreaStd*std(area);
end

if ~isfield(options,'erodeSize')
    erodeSize = round(sqrt(mean(area))/pi);
    disp(['separate fused nuclei erode size:' num2str(erodeSize)]);
else
    erodeSize = options.erodeSize;
end

if minSolidity > 0
    fusedCandidates = area./convexArea < minSolidity & area > minFusedArea;
else
    fusedCandidates = area > mean(area) + minAreaStd*std(area);
end
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1,sublist{:});

fusedMask = false(size(nuclearMask));
fusedMask(sublist) = 1;

nucmin = imerode(fusedMask,strel('disk',erodeSize));

outside = ~imdilate(fusedMask,strel('disk',1));
basin = imcomplement(bwdist(outside));
basin = imimposemin(basin, nucmin | outside);

L = watershed(basin);
newNuclearMask = nuclearMask & L > 1 | nuclearMask - fusedMask;

end