function [nucLevel,cytLevel,NCratio,bg] = readIntensityValues(img,masks,bgmask)
%for an image stack in a single channel, read out intensity values in
%nuclei, cytoplasm, and background

%replace zero-padding with NaNs before reading out values so padding is not
%included in averages; image must be cast to double first
if ~isa(img,'double')
    img = double(img);
end
img(img == 0) = NaN;

ncells = length(masks);
nz = size(img,3);

nucLevel = NaN(ncells,1);
cytLevel = NaN(ncells,1);
for cii = 1:ncells
    nucLevel(cii) = mean(img(masks(cii).nucmask),'omitnan');
    cytLevel(cii) = mean(img(masks(cii).cytmask),'omitnan');
end

bgmedians = NaN(nz,1);
for zi = 1:nz
    im = img(:,:,zi);
    bgmedians(zi) = median(im(bgmask(:,:,zi)),'omitnan');
end
bg = median(bgmedians,'omitnan');

NCratio = (nucLevel - bg)./(cytLevel - bg);

end



