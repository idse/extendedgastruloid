function zidxs = randomIndicesDistWeighted(upperleft, ims)
% find random indices for stitching overlapping images with a stitching
% weight from the image edge
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
M = imsize(1); N = imsize(2); %image size

totalSize = max(cat(1,upperleft{:})) + [M - 1, N - 1];

W = cell(1,1,4);
for ii = 1:4
    W{ii} = zeros(totalSize);
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
            
            
            upper = upperleft{ii,jj}(1);
            left = upperleft{ii,jj}(2);
            if ii == 1
                boverlap = M - (upperleft{ii+1,jj}(1) - upper); %bottom overlap
                leftweight = [ones(1,M-boverlap),linspace(1,0,boverlap)];
            elseif ii == m
                toverlap = M - (upper - upperleft{ii-1,jj}(1)); %top overlap
                leftweight = [linspace(0,1,toverlap),ones(1,M-toverlap)];
            else
                boverlap = M - (upperleft{ii+1,jj}(1) - upper);
                toverlap = M - (upper - upperleft{ii-1,jj}(1));
                leftweight = [linspace(0,1,toverlap),...
                    ones(1,M-toverlap-boverlap),linspace(1,0,boverlap)];
            end

            if jj == 1
                roverlap = N - (upperleft{ii,jj+1}(2) - left); %right overlap
                rightweight = [ones(1,N-roverlap),linspace(1,0,roverlap)];
            elseif jj == n
                loverlap = N - (left - upperleft{ii,jj-1}(2)); %left overlap
                rightweight = [linspace(0,1,loverlap),ones(1,N-loverlap)];
            else
                roverlap = N - (upperleft{ii,jj+1}(2) - left); %right overlap
                loverlap = N - (left - upperleft{ii,jj-1}(2)); %left overlap
                rightweight = [linspace(0,1,loverlap),...
                    ones(1,N-loverlap-roverlap),linspace(1,0,roverlap)];
            end
            
            W{Z(ii,jj)}(I,J) = leftweight'*rightweight;
        else
            disp(['skipping ' num2str([ii jj])]);
        end
    end
end
W = cell2mat(W);

[Wm,Wi] = max(W,[],3); %max weight and z index of max weight at each pixel
numinds = sum(W > 0,3); %number of nonzero weights at each pixel

%zidxs keeps track of the index to use at each pixel in the large image
zidxs = ones(size(Wi));
%there is only one option for most pixels; assign those indices here
mask = Wm == 1;
zidxs(mask) = Wi(mask);

%use (wieghted) random sampling for the other pixels
[row,col] = find(numinds > 1);
for ii = 1:length(row)
    ri = row(ii); ci = col(ii);
    p = squeeze(W(ri,ci,:));
    zidxs(ri,ci) = drawmnrnd(1,p); %draw from a multinomial distribution
end


end