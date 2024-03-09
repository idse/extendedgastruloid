function runPreprocessing(dataDir,stitchOptions,MIPoptions,zSliceOptions,...
    makeStitchedImages,makeMIPs,writeZslices)
%function to call separate preprocessing scripts given user inputs and
%metadata

% load(fullfile(dataDir,'meta.mat'),'meta');

%should there be any warning or error if more than one of these flags is
%set to true? in particular makeStitchedImages should generally be used
%alone

if makeMIPs
    inputdir = dataDir;
    outputdir = fullfile(dataDir, 'MIP');
    if ~exist(outputdir,'dir'), mkdir(outputdir); end
    
    batchMIP(inputdir, outputdir, MIPoptions);
end

if makeStitchedImages
    stitchImageMontages(dataDir, stitchOptions);
end

if writeZslices
    writeZstacks(dataDir,zSliceOptions)
end


end