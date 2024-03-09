function seg = ilastikRead(fname,binarize)

if ~exist('binarize','var')
    binarize = true;
end

seg = squeeze(h5read(fname,'/exported_data'));
if contains(fname, 'Simple Segmentation')
    if binarize
        seg = seg == 1;
    end
elseif contains(fname, 'Probabilities')
    %take first ilastik class and threshold if you have probabilities
    if binarize
        seg = squeeze(seg(1,:,:,:) > 0.5);
    else
        seg = squeeze(seg(1,:,:,:));
    end
elseif contains(fname, 'Object Predictions')
    %no extra step needed for object classification
else
    error('Unknown output type')
end

seg = permute(seg,[2,1,3]);

end
