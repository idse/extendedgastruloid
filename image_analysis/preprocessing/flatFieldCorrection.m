function newims = flatFieldCorrection(ims,G,D)
Dmean = mean(D,'all');

newims = cell(size(ims));
for ii = 1:size(ims,1)
    for jj = 1:size(ims,2)
        newims{ii,jj} = uint16((double(ims{ii,jj}) - D).*G + Dmean);
    end
end

end
