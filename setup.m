%-----------------------------------------
% set up path
%-----------------------------------------

% entries of matlab search path
pathentries = regexp(path, pathsep, 'split');

% find ones containing ImSAnE and remove from path
disp('removing old entries from path');
for i = 1:length(pathentries)
    if ~isempty(regexp(pathentries{i}, 'HeemskerkLab','once'))...
       || ~isempty(regexp(pathentries{i}, 'heemskerklab','once'))
        rmpath(pathentries{i});
    end
end

%%
% path of current script is new ImSAnE directory
[repopath,~,ext] = fileparts(matlab.desktop.editor.getActiveFilename);
disp(['repo path: ' repopath]);
addpath(genpath(repopath));
disp('added directory containing setup.m and subdirectories to path');

% storing the path in the startup directory
upath = userpath;
savepath(fullfile(upath, 'pathdef.m'));
