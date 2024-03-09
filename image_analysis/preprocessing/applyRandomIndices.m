function stitched = applyRandomIndices(upperleft, ims, zidxs)
% apply random indices in overlap regions to stitch a large image
%
% upperleft:    cell array of positions of upperleft corners with
%               upperleft corner of upperleft image being (1,1)
%               as produced by registerImageGrid(..)
% imgs:         cell array of images 
%               with rows and cols corresponding to grid

UL = cat(1,upperleft{:});
dUL = 1 - min(UL);
for ii = 1:numel(upperleft)
    if ~isempty(upperleft{ii})
        upperleft{ii} = upperleft{ii} + dUL;
    end
end

% not all positions on grid have to contain an image
% so find the first actual image on the grid
emptycells = cellfun(@isempty,ims,'UniformOutput',false);
emptycells = cat(1,emptycells{:});
cidx = find(~emptycells,1,'first');
imsize = size(ims{cidx},[1,2]);
imclass = class(ims{cidx}); %data type of the images (usually uint16)
M = imsize(1); N = imsize(2); %image size
totalSize = max(cat(1,upperleft{:})) + [M - 1, N - 1];

IS = cell(1,1,4);
for ii = 1:4
    IS{ii} = zeros(totalSize,imclass);
end

m = size(ims,1);
n = size(ims,2);
Z = zeros(m,n);
starts = [1,1;1,2;2,1;2,2];
for idx = 1:4
    for ii = starts(idx,1):2:m
        for jj = starts(idx,2):2:n
            Z(ii,jj) = idx;
        end
    end
end

for ii = 1:m
    for jj = 1:n
        if ~isempty(ims{ii,jj})            
            I = upperleft{ii,jj}(1):upperleft{ii,jj}(1)+M-1;
            J = upperleft{ii,jj}(2):upperleft{ii,jj}(2)+N-1;
            
            IS{Z(ii,jj)}(I,J) = ims{ii,jj};
        else
            disp(['skipping ' num2str([ii jj])]);
        end
    end
end

stitched = zeros(totalSize,imclass);
for ii = 1:4
    mask = zidxs == ii;
    stitched(mask) = IS{ii}(mask);
end





end