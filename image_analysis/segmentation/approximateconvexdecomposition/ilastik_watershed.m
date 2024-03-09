function nucseg = ilastik_watershed(probs,opts)
%function to separate touching nuclei in an Ilastik probability map using
%hysteresis thresholding and watershed. The core threshold is used to find
%objects, i.e. each region of connected pixels exceeding the core (higher)
%threshold is the marker of a single object. These regions are used as
%seeds for a seeded watershed of the complement of the probability mask.
%finally, background pixels are removed with the final threshold

if numel(size(probs)) > 2
    probs = squeeze(probs(1,:,:));
end

%parse inputs
if ~exist('opts','var')
    opts = struct();
end
if ~isfield(opts,'core' )
    %detection probability threshold
    opts.core = 0.8;
end
if ~isfield(opts,'final')
    %segmentation probability threshold
    opts.final = 0.5;
end
if ~isfield(opts,'minSize')
    %minimum size of an acceptable obect seed
    opts.minSize = 10;
end
if ~isfield(opts,'maxSize')
    %maximum size of an acceptable obect seed*****
    opts.maxSize = 100000;
end

core = opts.core; final = opts.final;
minSize = opts.minSize; maxSize = opts.maxSize;
%find markers for individual objects with the core threshold
seeds = probs > core;
%remove very large and very small objects from the detection mask
CC = bwconncomp(seeds);
sizes = cellfun(@numel, CC.PixelIdxList);
idxs = CC.PixelIdxList(sizes < minSize | sizes > maxSize);
idxs = cell2mat(idxs(:));
seeds(idxs) = 0;
%perform marker-controlled watershed of complement of probability mask
V = imimposemin(imcomplement(probs), seeds);
nucseg = watershed(V);

%determine pixel-wise segmentation mask via thresholding
mask = probs > final;
%remove background pixels from the watershed result
nucseg(~mask) = 0;
%make the watershed result a binary segmentation instead of a label matrix
nucseg = nucseg > 0;

end