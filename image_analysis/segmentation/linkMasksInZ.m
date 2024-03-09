function [masks, cellData, chain] = linkMasksInZ(allmasks, meta, opts)
% create new cellData merging info from different zslices

%options
if ~isfield(opts,'minZsize')
    opts.minZsize = 1;
end

if ~isfield(opts,'suppressOutput')
    suppressOutput = false;
else
    suppressOutput = opts.suppressOutput;
end

% parfor doesn't work with struct arrays, forcing us to use cell arrays
if iscell(allmasks)
    disp('converting cell array to struct array')
    allmasks = cellofstructs2structarray(allmasks);
end

%IoU = intersection over union threshold (sort of)
IoU = opts.IoU;
%convert microns to pixels or frames for other options
maxZslices = max([floor(opts.maxZsize/meta.zres),1]);
minZslices = floor(opts.minZsize/meta.zres);
maxCentroidDistance = opts.maxCentroidDistance/meta.xres;

% match nuclei between frames
nz = length(allmasks);
% determine first non-empty slice; nei = non-empty index
nei = find(cellfun(@(x) ~isempty(x),{allmasks.nucmask}),1,'first');

if ~suppressOutput
    fprintf('Max z slices = %d slices = %g um\n',...
        maxZslices,maxZslices*meta.zres)
    fprintf('Min z slices = %d slices = %g um\n',...
        minZslices,minZslices*meta.zres)
    fprintf('max centroid distance = %g um = %g pixels\n',...
        maxCentroidDistance*meta.xres,maxCentroidDistance)
    fprintf('First non-empty z index = %d\n',nei)
end

%number of cells in the first frame; use this to keep track of the number
%of uniquely identified nuclei across slices
cn = length(allmasks(nei).nucmask); %cn = cell number
%build chains of linked nuclear masks
chain = cell(cn,1);
for jj = 1:cn
    %chain{jj} = [z slice, cell index, nuclear area]
    nucArea = numel(allmasks(nei).nucmask{jj});
    chain{jj} = [nei,jj,nucArea];
end
%for each slice, keep track of the chain to which each cell is assigned;
%for the first z slice, each cell initializes its own chain
cidxs = cell(nz,1);
cidxs{nei} = (1:cn)';

%do the assignment as a linear assignment problem between adjacent z-slices
%with costs based on how closely objects overlap
for ii = nei:nz-1
    XY1 = allmasks(ii).XY;
    XY2 = allmasks(ii+1).XY;
    %number of cells in each frame
    nc1 = size(XY1,1); nc2 = size(XY2,1);
    %don't try to link if one or both frames have no cells
    if (nc1 > 0) && (nc2 > 0)
        % impose minimal nubmer of cells (crash if ncells < k)?
        k = min([nc1,nc2,3]);
        A11 = Inf(nc1,nc2);
        
        %only calculate overlap for a few possible candidates for each nucleus
        %but do it both ways (neighbors in slice ii+1 of a cell in slice ii and
        %neighbors in slice ii of a cell in slice ii+1)
        [Idx, D] = knnsearch(XY1,XY2,'K',k);
        for jj = 1:nc2
            list2 = allmasks(ii+1).nucmask{jj};
            for ki = 1:k
                ll = Idx(jj,ki);
                if D(jj,ki) < maxCentroidDistance
                    list1 = allmasks(ii).nucmask{ll};
                    A11(ll,jj) =...
                        min(numel(list1),numel(list2))/numel(intersect(list1,list2));
                end
            end
        end

        [Idx, D] = knnsearch(XY2,XY1,'K',k);
        for jj = 1:nc1
            list1 = allmasks(ii).nucmask{jj};
            for ki = 1:k
                ll = Idx(jj,ki);
                if D(jj,ki) < maxCentroidDistance
                    list2 = allmasks(ii+1).nucmask{ll};
                    A11(jj,ll) =...
                        min(numel(list1),numel(list2))/numel(intersect(list1,list2));
                end
            end
        end

        %make alternative cost matrices to reject links using the IoU threshold
        %Inf square matrix with alternative costs on the diagonal
        A12 = Inf(nc1,nc1);
        A12(eye(nc1)==1) = 1/IoU;
        A21 = Inf(nc2,nc2);
        A21(eye(nc2)==1) = 1/IoU;

        %final cost matrix
        CM = [...
            A11, A12;...
            A21, A11'];
        %do the optimization
        [CM_indices, ~] = lapjv(CM);
        %parse resulting assignments
        target_indices = CM_indices(1:nc1);
        source_indices = CM_indices((nc1+1):(nc1+nc2)) - nc2;
        %link rejection is denoted by index of -1
        target_indices(target_indices > nc2) = -1;
        source_indices(source_indices < 1) = -1;

        %cells in frame ii+1 that were assigned to a cell in frame ii get added
        %to the chain that that cell is a part of
        C = NaN(nc2,1);
        for jj = 1:nc1
            tidx = target_indices(jj);
            cidx = cidxs{ii}(jj);
            if tidx > 0
                nucArea = numel(allmasks(ii+1).nucmask{tidx});
                chain{cidx} = [chain{cidx};ii+1,tidx,nucArea];
                C(tidx) = cidx;
            end
        end

        %nuclei in frame ii+1 with no link to one in frame ii get new chains
        %instead of being assigned to existing ones
        newinds = find(source_indices < 0);
        for jj = 1:length(newinds)
            cn = cn + 1;
            tidx = newinds(jj);
            nucArea = numel(allmasks(ii+1).nucmask{tidx});
            chain{cn} = [ii+1,tidx,nucArea];
            C(tidx) = cn;
        end
        cidxs{ii+1} = C;
        
    else
        %if there are no cells in frame ii or ii + 1, every cell in frame
        %ii+1 gets a new chain
        cidxs{ii+1} = (cn+1):(cn+nc2);
        for jj = 1:nc2
            cn = cn + 1;
            nucArea = numel(allmasks(ii+1).nucmask{jj});
            chain{cn} = [ii+1,jj,nucArea];
        end
    end
end

%deal with chains with more than maxZslices components
I = find(cellfun(@(x) size(x,1), chain) > maxZslices);
while ~isempty(I)
    for ii = 1:length(I)
        areas = chain{I(ii)}(:,3);
        areaDiffs = areas(2:end) - areas(1:end-1);
        secondDiffs = zeros(size(areas));
        secondDiffs(2:end-1) = areaDiffs(2:end) - areaDiffs(1:end-1);
        filler = min(secondDiffs)-1;
        secondDiffs([1,end]) = filler;
        [~,Idx] = max(secondDiffs);
        c1 = chain{I(ii)}(1:Idx,:);
        c2 = chain{I(ii)}(Idx+1:end,:);
        chain{I(ii)} = c1;
        cn = cn + 1;
        chain{cn} = c2;
    end
    I = find(cellfun(@(x) size(x,1), chain) > maxZslices);
end
%filter out chains with less than minZslices components
chain = chain(cellfun(@(x) size(x,1), chain) >= minZslices);

% reduce cellData to have info for only one z slice for each cell
%make each nucleus have one set of pixels, defined in the slice in which
%its average intensity is highest
nCells = length(chain);
masks(nCells) = struct;
cellData = struct;
m = meta.ySize; n = meta.xSize;

%initialize fields in cellData
fields = fieldnames(allmasks(nei).nucstats);
cellData.XY = NaN(nCells,2);
cellData.nucArea = NaN(nCells,1);
cellData.nucVolume = NaN(nCells,1);
cellData.nucZ = NaN(nCells,1); %z slice of largest nuclear area
cellData.Z = NaN(nCells,1); %centroid in z
nfields = length(fields);
for fi = 1:nfields
    cellData.(fields{fi}) = NaN(nCells,1);
end

%for each chain, find the slice for which the nucleus is brightest in the 
%nuclear channel and use data from that slice for the cell
for jj = 1:nCells
    %cellinfo(ii,:) = [z slice, cell index, nuclear area]
    cellinfo = chain{jj};
    
    %combine the slice-wise pixel index lists into a list of (linear)
    %indices for all the voxels occupied by the cell in 3D
    pixelindices = cell(size(cellinfo,1),1);
    for zi = 1:size(cellinfo,1)
        z = cellinfo(zi,1); % z slice
        ci = cellinfo(zi,2); %cell index
        pixelindices{zi} = allmasks(z).nucmask{ci} + (z - 1)*m*n;
    end
    %find cell centroid in z
    npixels = cellfun(@numel,pixelindices);
    cellData.Z(jj) = sum(npixels.*cellinfo(:,1))/sum(npixels);
    cellData.nucVolume(jj) = sum(npixels);
    
    masks(jj).nucmask = cell2mat(pixelindices);
    
    %only keep the cytmask in the z slice with the largest nuclear area
    [~, I] = max(cellinfo(:,3)); %row of cellinfo with highest nuclear level
    zi = cellinfo(I,1); %z slice of largest nuclear size
    ci = cellinfo(I,2); %cell index in that z slice
    if ~isempty(allmasks(zi).cytmask)
        masks(jj).cytmask = allmasks(zi).cytmask{ci} + (zi - 1)*m*n;
    else
        masks(jj).cytmask = [];
    end
    cellData.nucZ(jj) = zi;
    
    cellData.nucArea(jj) = cellinfo(I,3);
    for fi = 1:nfields
        cellData.(fields{fi})(jj) = allmasks(zi).nucstats(ci).(fields{fi});
    end
    cellData.XY(jj,:) = allmasks(zi).XY(ci,:);
end





end