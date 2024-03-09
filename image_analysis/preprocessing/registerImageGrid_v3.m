function [upperleft,links] = registerImageGrid_v3(imgs, pixelOverlap,maxshift)
    % register a grid of overlapping images
    % 
    % [upperleft,links] = registerImageGrid(imgs, pixelOverlap)
    %
    %
    % upperleft:    cell array of positions of upperleft corners with
    %               upperleft corner of upperleft image being (1,1)
    % links:        a matrix containing correlation links between images
    %               empty positions will correspond to cells themselves
    %
    % imgs:         cell array of images 
    %               with rows and cols corresponding to grid
    % pixelOverlap: width of overlapping strip in pixels
    %               if left empty, upperleft for a square grid of images
    %               will be returned
    
    
    
    % ---------------------
    % Idse Heemskerk, 2016
    % v2, Seth Teague 2020: solve optimization problem to use redundant
    % shifts to improve robustness
    % updated Seth Teague 2022: allow non-square images on a non-square
    % grid
    % ---------------------
    
%maximum allowable shift in x or y so large shifts detected due to noise in
%overlap regions that contain no cells are not used
if ~exist('maxshift','var')
    maxshift = 20;
end

% find the first actual image on the grid
emptycells = cellfun(@isempty,imgs,'UniformOutput',false);
emptycells = cat(1,emptycells{:});
imsize = size(imgs{find(~emptycells,1,'first')},[1,2]);
M = imsize(1); N = imsize(2);

if numel(pixelOverlap) == 1
    pixelOverlapY = pixelOverlap;
    pixelOverlapX = pixelOverlap;
elseif numel(pixelOverlap) == 2
    pixelOverlapY = pixelOverlap(1); pixelOverlapX = pixelOverlap(2);
end

% make upperleft for square grid of images (no stitching)
grid = size(imgs,[1,2]);
gridUpperleft = cell(size(imgs));
for i = 1:size(imgs, 1)
    for j = 1:size(imgs, 2)
        gridUpperleft{i,j} = [1+(i-1)*(M + 50), 1+(j-1)*(N + 50)];
    end
end

m = grid(1); n = grid(2);

%shifts between image (i,j) and image (i+1,j)
len = (m - 1)*n;
b1_x = zeros(len,1); b1_y = zeros(len,1); A1 = zeros(len,m*n);
for ii = 1:m-1
    for jj = 1:n
        if ~isempty(imgs{ii,jj}) && ~isempty(imgs{ii+1,jj})
            img = imgs{ii,jj}(end-pixelOverlapY+1:end,:);
            below = imgs{ii+1,jj}(1:pixelOverlapY,:);
            [shifty, shiftx, ~, ~] = xcorr2fft(img, below);
            
            if abs(shiftx) > maxshift || abs(shifty) > maxshift
                shiftx = 0; shifty = 0;
            end
        else
            shifty = 0; shiftx = 0;
        end
        idx = (m-1)*(jj-1) + ii;
        b1_y(idx) = M - pixelOverlapY - shifty;
        b1_x(idx) = -shiftx;
        
        idx1 = m*(jj-1) + ii;
        idx2 = m*(jj-1) + ii + 1;
        A1(idx,idx1) = -1; A1(idx,idx2) = 1;
    end
end

%shifts between image (i,j) and image (i,j+1)
len = (n - 1)*m;
b2_x = zeros(len,1); b2_y = zeros(len,1); A2 = zeros(len,m*n);
for ii = 1:m
    for jj = 1:n-1
        if ~isempty(imgs{ii,jj}) && ~isempty(imgs{ii,jj+1})
            img = imgs{ii,jj}(:,end-pixelOverlapX+1:end);
            right = imgs{ii,jj+1}(:,1:pixelOverlapX);
            [shifty, shiftx, ~, ~] = xcorr2fft(img, right);
            
            if abs(shiftx) > maxshift || abs(shifty) > maxshift
                shiftx = 0; shifty = 0;
            end
            
        else
            shifty = 0; shiftx = 0;
        end
        idx = (n-1)*(ii-1) + jj;
        b2_y(idx) = -shifty;
        b2_x(idx) = N - pixelOverlapX - shiftx;
        
        idx1 = m*(jj-1) + ii;
        idx2 = m*jj + ii;
        A2(idx,idx1) = -1; A2(idx,idx2) = 1;
    end
end

A = [A1; A2]; bx = [b1_x; b2_x]; by = [b1_y; b2_y];
%take the matrix pseudoinverse for least squares solution
pA = pinv(A);
xposes = pA*bx;
yposes = pA*by;
%define the position of the top left image to be (0,0)
xposes = xposes - xposes(1);
yposes = yposes - yposes(1);

xpositions = round(reshape(xposes,[m,n]));
ypositions = round(reshape(yposes,[m,n]));

upperleft = cell(grid);
for ii = 1:grid(1)
    for jj = 1:grid(2)
        upperleft{ii,jj} = [ypositions(ii,jj), xpositions(ii,jj)];
    end
end
