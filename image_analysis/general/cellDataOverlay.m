function [colormap, colidx] = cellDataOverlay(im, XY, celldata, tol)

    [N,bins] = histcounts(celldata,20);
    cumdist = cumsum(N/sum(N));
    maxval = bins(find(cumdist > 1-tol, 1, 'first')); 
    minval = bins(find(cumdist > tol, 1, 'first')); 
    
    celldata(celldata > maxval) = maxval;
    celldata(celldata < minval) = minval;
    
    celldata = (celldata - minval)/(maxval - minval); 

    colidx = floor(celldata*255 + 1);
    colormap = jet(256);

    if size(im,3) == 1
        im = mat2gray(im);
        A = imadjust(im,stretchlim(im,0.001));
    else
        A = im;
    end
    
    imshow(A);
    hold on
    scatter(XY(~isnan(celldata),1), XY(~isnan(celldata),2),20,colormap(colidx(~isnan(celldata)),:),'filled')
    hold off
end