function [stats, trueNucInfo, nucInfo] =...
    detection_accuracy(mask, true_mask, opts)

%find objects in ground truth and test masks, filter out very small (noise)
%objects, define a cutoff distance based on average nuclear radius in
%ground truth mask, construct a distance-based cost matrix for linking
%objects between one another, assigning cost of Inf for nuclei farther
%apart than the desired distance
%to filter out objects that are in Ilastik base mask but not in ground
%truth (missegmentation at the pixel level)-> mask = mask & true_mask?

if ~isfield(opts,'minArea')
    minArea = 50;
else
    minArea = opts.minArea;
end

%if the mask has extraneous color channels, remove them (all relevant info
%is in the red channel)
if ndims(true_mask) == 3
    true_mask = squeeze(true_mask(:,:,1));
end

%remove objects that are in the calculated segmentation but not the gold
%truth segmentation (pixel-level classification errors, we are interested
%in object-level)
mask = mask & true_mask;

%initialize arrays to hold cell data
area = [];
idxlists = {};
centroids = [];
diameters = [];

%get data on nuclei from each object in the true mask
vals = unique(true_mask);
vals = vals(vals > 0);
for vi = 1:length(vals)
    CC = bwconncomp(true_mask == vals(vi), 4);
    props = regionprops(CC,'Area', 'Centroid','EquivDiameter');
    areas = {props.Area};
    areas = cell2mat(areas(:));
    diameter = {props.EquivDiameter};
    diameter = cell2mat(diameter(:));
    centroid = {props.Centroid};
    centroid = cell2mat(centroid(:));
    
    area = [area; areas];
    idxlists = [idxlists, CC.PixelIdxList];
    diameters = [diameters; diameter];
    centroids = [centroids; centroid];
end

kept = area > minArea;
tempseg = zeros(size(mask));
allidxs = idxlists(kept); allidxs = cell2mat(allidxs(:));
tempseg(allidxs) = 1;

cutoff = 0.5*mean(diameters);

%find the number of objects that are inside clumps in the true mask
CC = bwconncomp(tempseg);
%if a connected component in the true mask has only one value shared
%between all pixels, then it is not a clump; if we remove all these single,
%unclumped objects, we know how many objects are within clumps
values = cellfun(@(x) numel(unique(true_mask(x))), CC.PixelIdxList);
nclumps = sum(values > 1);
nclumped = sum(kept) - sum(values == 1);
nsplits = nclumped - nclumps;

trueNucInfo = struct(...
    'area',         area(kept),...
    'idxs',         {idxlists(kept)},...
    'centroid',     centroids(kept,:),...
    'radius',       diameters(kept)/2);

CC = bwconncomp(mask > 0, 4);
props = regionprops(CC,'Area', 'Centroid','EquivDiameter');
areas = {props.Area};
areas = cell2mat(areas(:));
diameter = {props.EquivDiameter};
diameter = cell2mat(diameter(:));
centroid = {props.Centroid};
centroid = cell2mat(centroid(:));

kept = areas > minArea;

nucInfo = struct(...
    'area',         areas(kept),...
    'idxs',         {CC.PixelIdxList(kept)},...
    'centroid',     centroid(kept,:),...
    'radius',       diameter(kept)/2);

%set up linear assignment problem to match centroids
source = trueNucInfo.centroid;
target = nucInfo.centroid;
n_source_points = size(source, 1);
n_target_points = size(target, 1);
A1 = NaN(n_source_points, n_target_points);
%build distance-based cost matrix, with distances higher than cutoff set to
%Inf
for i = 1 : n_source_points
    % Pick one source point
    current_point = source(i, :); %xy value of current point

    % Compute square distance to all target points
    diff_coords = target - repmat(current_point, n_target_points, 1);
    square_dist = sum(diff_coords.^2, 2);
    square_dist(square_dist > cutoff^2) = Inf;
    
    % Store them
    A1(i, :) = square_dist;
end


linking_costs = A1(~isinf(A1)); %extract the finite costs
thresh = 1.1*max(linking_costs,[],'all'); %use max value
A2 = thresh*eye(n_source_points);
A2(A2 == 0) = Inf;

A3 = thresh*eye(n_target_points);
A3(A3 == 0) = Inf;

CM = [A1, A2; A3, transpose(A1)];

% error('whats going on')
[CM_indices, ~] = lapjv(CM);

target_indices = CM_indices(1:n_source_points);
source_indices = CM_indices((n_source_points+1):(n_source_points+n_target_points)) - n_target_points;
target_indices(target_indices > n_target_points) = -1;
source_indices(source_indices < 1) = -1;

%true positives are objects that are matched up in both
TP = find(target_indices > 0);
%false positives are objects in mask (target) that are not in true_mask (source)
FP = find(source_indices == -1);
%false negatives are objects in true_mask (source) that are not in mask (target)
FN = find(target_indices == -1);

%find dice similarity coefficient, jaccard index for true positives
%JI(A,B) = |intersection(A,B)|/|union(A,B)| -> ranges from 0 to 1
%DSC(A,B) = 2*|intersection(A,B)|/(|A| + |B|) -> ranges from 0 to 1
%indices of pixels in an object
JI = zeros(length(TP),1);
DSC = zeros(length(TP),1);
for ti = 1:length(TP)
    sidx = TP(ti); %source (true) index
    tidx = target_indices(sidx); %target (calculated) index
    A = trueNucInfo.idxs{sidx}; %pixels in source
    B = nucInfo.idxs{tidx}; %pixels in target
    JI(ti) = numel(intersect(A,B))/numel(union(A,B));
    DSC(ti) = 2*numel(intersect(A,B))/(numel(A) + numel(B));
end

%new calculation for true positives: (it would be better to actually
%identify which objects are inside clumps instead of doing it this way)
%nTP = numel(TP); nFP = numel(FP); nFN = numel(FN);
%we want to evaluate how many splits were made correctly or incorrectly;
%each false negative object detection indicates a split that should have
%been made and was not; each false positive indicates a split that was made
%erroneously
nFP = numel(FP); nFN = numel(FN); nTP = nsplits - nFN;

%calculate precision, recall, accuracy, F-measure
%Precision = TP/(TP + FP)
%Recall = TP/(TP + FN)
%Accuracy = TP/(TP + FN + FP)
%F-measure = 2*(Precision*Recall)/(Precision + Recall)

Precision = nTP/(nTP + nFP);
Recall = nTP/(nTP + nFN);
Accuracy = nTP/(nTP + nFP + nFN);
Fmeasure = 2*Precision*Recall/(Precision + Recall);

%store true positives, false positives, false negatives, dice similarity
%coefficients and jaccard indices
stats = struct(...
    'TP',   TP',...
    'FP',   FP',...
    'FN',   FN',...
    'DSC',  DSC,...
    'JI',   JI,...
    'Precision',    Precision,...
    'Recall',       Recall,...
    'Accuracy',     Accuracy,...
    'Fmeasure',     Fmeasure,...
    'centroid',     nucInfo.centroid,...
    'trueCentroid', trueNucInfo.centroid,...
    'nclumped',     nclumped);

end



