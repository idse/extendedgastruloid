function saveTiff(im, fname, mode)

    if isa(im,'uint16')
        nbits = 16;
    else
        error('add other cases');
    end

    t = Tiff(fname,mode);
    
    numrows = size(im,1);
    numcols = size(im,2);
    
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.ImageLength = numrows;
    tagstruct.BitsPerSample = nbits;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.ImageWidth = numcols;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    
    t.setTag(tagstruct);
    t.write(im(:,:,1));
    for i=2:size(im,3)
        t.writeDirectory();
        t.setTag(tagstruct);
        t.write(im(:,:,i));
    end
    t.close();   
end