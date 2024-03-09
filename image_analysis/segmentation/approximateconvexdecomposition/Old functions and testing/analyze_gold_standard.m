%% load data
clear; close all; clc;
dataDir = 'G:\My Drive\Research\Heemskerk lab\Data\convex_decomposition_gold_standard';
load(fullfile(dataDir,'gold_data.mat'));
diagnostics = gold_data.diagnostics;

%%
close all
img = gold_data.img';
seg = gold_data.seg';
newseg = seg;
CC = bwconncomp(seg);

badpixels = CC.PixelIdxList([diagnostics.bad] == 1);
badpixels = cell2mat(badpixels(:));
newseg(badpixels) = 0;

badidxs = find([diagnostics.bad] == 1);
imshow(cat(3, img, 0.7*img + 0.3*seg, 0.7*img + 0.3*seg))

for bi = 1:length(badidxs)
    bounds = diagnostics(badidxs(bi)).boundingbox;
    yl = bounds(1,:);
    xl = bounds(2,:);
    line([yl(1);yl(1);yl(2);yl(2);yl(1)],...
        [xl(1);xl(2);xl(2);xl(1);xl(1)],'Color','r')
    title(num2str(bi))
    pause(0.125)
end

dividing = CC.PixelIdxList([diagnostics.dividing] == 1);
dividing = cell2mat(dividing(:));
newseg(dividing) = 0;
% imshowpair(seg,newseg)

dividxs = find([diagnostics.dividing] == 1);
for bi = 1:length(dividxs)
    bounds = diagnostics(dividxs(bi)).boundingbox;
    yl = bounds(1,:);
    xl = bounds(2,:);
    line([yl(1);yl(1);yl(2);yl(2);yl(1)],...
        [xl(1);xl(2);xl(2);xl(1);xl(1)],'Color','g')
    title(num2str(bi))
    pause(0.125)
end

% manidxs = find([diagnostics.manual] == 1);
% for bi = 1:length(dividxs)
%     bounds = diagnostics(dividxs(bi)).boundingbox;
%     yl = bounds(1,:);
%     xl = bounds(2,:);
%     line([yl(1);yl(1);yl(2);yl(2);yl(1)],...
%         [xl(1);xl(2);xl(2);xl(1);xl(1)],'Color','g')
%     title(num2str(bi))
%     pause(0.125)
% end

% imshow(double(cat(3, seg, newseg, newseg)))

%%
close all
clc
splits = find([diagnostics.split] == 1);
notsplit = find([diagnostics.split] == 0);
manuals = find([diagnostics.manual] == 1);
notmanual = find([diagnostics.manual] == 0);
dividing = find([diagnostics.dividing] == 1);
bads = find([diagnostics.bad] == 1);

%good indices do not include bad or dividing cells
gidxs = (1:length(diagnostics))';
gidxs = gidxs(~ismember(gidxs,bads));
gidxs = gidxs(~ismember(gidxs,dividing));

%discard bad and dividing cells
splits = splits(ismember(splits,gidxs));
notsplit = notsplit(ismember(notsplit,gidxs));
notmanual = notmanual(ismember(notmanual,gidxs));

%false negatives are cells that I identified as needing to be split but
%which were not split by the algorithm
false_negatives = manuals(~ismember(manuals, splits));
%false positives are cells I said should not be split but which were anyway
false_positives = splits(~ismember(splits, manuals));
%true positives are cells which should have been split and were
true_positives = manuals(ismember(manuals,splits));
%true negatives are cells which should not have been split and were not
true_negatives = notmanual(ismember(notmanual,notsplit));

disp(strcat("Number of false negatives = ", num2str(numel(false_negatives))))
disp(strcat("Number of false positives = ", num2str(numel(false_positives))))
disp(strcat("Number of true negatives = ", num2str(numel(true_negatives))))
disp(strcat("Number of true positives = ", num2str(numel(true_positives))))


%Plot absolute concavity vs major axis length for merged and non-merged
%cases
% subplot(1,2,1)
curvcat = find(([diagnostics.curve1] < 80).*([diagnostics.curve2] < 80));
%cells that divide and meet curvature requirement
cat1 = intersect(manuals,curvcat);
%cells that divide and do not meet curvature requirement
cat2 = manuals(~ismember(manuals,curvcat));
%cells that do not divide and meet curvature requirement
cat3 = intersect(notmanual,curvcat);
%cells that do not divide and do not meet curvature requirement
cat4 = notmanual(~ismember(notmanual,curvcat));

scatter([diagnostics(cat1).concave1],...
    [diagnostics(cat1).concave2],36,'c','filled')
hold on
scatter([diagnostics(cat2).concave1],...
    [diagnostics(cat2).concave2],36,'b','filled')
scatter([diagnostics(cat3).concave1],...
    [diagnostics(cat3).concave2],36,'m','filled')
scatter([diagnostics(cat4).concave1],...
    [diagnostics(cat4).concave2],36,'r')
title('Relative vs absolute concavity for manually annotated nuclei')
xlabel('Absolute concavity')
xline(5);
yline(1);
ylabel('Concavity to cut length ratio')
legend('split, curve', 'split no curve', 'no split curve',...
    'no split no curve','tau1','tau2')
% subplot(1,2,2)
% title('Should not be split')
% xlabel('Absolute concavity')
% ylabel('Major Axis Length')
set(gcf,'WindowState','maximized')
set(gca,'yscale','log')

%% Save images of false negatives
close all
saveDir = fullfile(dataDir,"false_negatives");
pattern = 'example_t%04d.tif';
list = dir(saveDir);
imshow(cat(3,img,0.5*img + 0.5*seg,img))
for ti = 1:length(false_negatives)
    fprintf('.')
    idx = false_negatives(ti);
    stuff = diagnostics(idx);
    xl = stuff.boundingbox(1,:);
    yl = stuff.boundingbox(2,:);
    height = xl(2) - xl(1);
    xl = xl + [-height/2, height/2];
    width = yl(2) - yl(1);
    yl = yl + [-width/2, width/2];
    xlim(xl)
    ylim(yl)
    title(strcat("Absolute concavity = ", num2str(stuff.concave1),...
        ", Relative concavity = ", num2str(stuff.concave2)))
    xlabel(strcat("Angle 1 = ", num2str(stuff.curve1), ", Angle 2 = ",...
        num2str(stuff.curve2)))
    set(gcf,'WindowState','maximized')
    frame = getframe(gcf);
    writename = fullfile(saveDir,sprintf(pattern,ti));
    imwrite(frame.cdata,writename)
%     pause
end
fprintf('\n')

%% Plot absolute concavity vs relative concavity for merged and non-merged cases
close all
scatter([diagnostics(manuals).concave1],...
    [diagnostics(manuals).concave2],36,'c','filled')
hold on
scatter([diagnostics(notmanual).concave1],...
    [diagnostics(notmanual).concave2],36,'m')
title('Realative vs absolute concavity for manually annotated nuclei')
xlabel('Absolute concavity')
ylabel('Relative concavity')
legend('Should be split','Should not be split')
set(gcf,'WindowState','maximized')
set(gca,'yscale','log')

%% Examine objects that should be split
close all
imshow(img)
for ti = 1:length(manuals)
    idx = manuals(ti);
    xl = diagnostics(idx).boundingbox(1,:); width = (xl(2) - xl(1))/2;
    yl = diagnostics(idx).boundingbox(2,:); height = (yl(2) - yl(1))/2;
    xl = [xl(1) - width, xl(2) + width];
    yl = [yl(1) - height, yl(2) + height];
    xlim(xl);
    ylim(yl)
    c1 = diagnostics(idx).concave1;
    c2 = diagnostics(idx).concave2;
    c3 = diagnostics(idx).curve1;
    c4 = diagnostics(idx).curve2;
    title(strcat("Abs concavity = ", num2str(c1),...
        "; Rel concavity = ", num2str(c2)))
    xlabel(strcat("Curve 1 = ", num2str(c3), "; Curve 2 = ", num2str(c4)))
    disp(ti)
    pause
end

%% Make ROC curves
close all
manuals = [diagnostics.manual] == 1;
dividing = [diagnostics.dividing] == 1;
bads = [diagnostics.bad] == 1;
% discard = dividing | bads;
discard = bads;
keep = ~discard;
tau1s = linspace(1,25,25);
tau2s = linspace(0.25,2,25);
concave1 = [diagnostics.concave1];
concave2 = [diagnostics.concave2];
curvcat = ([diagnostics.curve1] < 80).*([diagnostics.curve2] < 80);

figure
%ROC for tau2 with tau1 fixed at 5
tps = zeros(25,5);
fps = zeros(25,5);
tau1s = [0.5, 1, 5, 10, 15];
leg = {'tau1 = 0.5', 'tau1 = 1','tau1 = 5','tau1 = 10','tau1 = 15'};
for ji = 1:5
    tau1 = tau1s(ji);
    for ti = 1:25
        tau2 = tau2s(ti);
        split = (concave1 > tau1).*(concave2 > tau2 | curvcat);
        split = split.*keep;
        npos = sum(split);
        %true positives
        tp = split.*manuals;
        tps(ti,ji) = sum(tp);
        %false positives
        fp = split.*(~manuals);
%         fps(ti,ji) = npos - tps(ti,ji);
        fps(ti,ji) = sum(fp);
    end
end
tps = 100*tps/sum(manuals);
fps = 100*fps/sum(manuals);
plot(fps,tps)
xlabel("False positive rate")
ylabel("True positive rate")
legend(leg,'Location','southeast')
xlim([0,150])
ylim([60, 105])
title("ROC curve while varying relative concavity threshold from 0.25 to 2")

%%
close all
xl = [1343 1614];
yl = [3146 3417];
testseg = seg(yl(1):yl(2),xl(1):xl(2));
B = bwboundaries(testseg);
b = B{1};

Vertices = b;
Lines = [(1:length(Vertices)-1)', (2:length(Vertices))'];
tic
k=LineCurvature2D(Vertices,Lines);
N=LineNormals2D(Vertices,Lines);
toc

figure,  hold on;
k=k*10;
plot([Vertices(:,1) Vertices(:,1)+k.*N(:,1)]',[Vertices(:,2) Vertices(:,2)+k.*N(:,2)]','g');
plot([Vertices(Lines(:,1),1) Vertices(Lines(:,2),1)]',[Vertices(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
plot(Vertices(:,1),Vertices(:,2),'r.');
axis equal;
pause

figure
line(b(:,1),b(:,2))
npoints = size(Vertices,1);
vecs = zeros(size(N));
scaled = N.*k;
weights = [1; 2; 3; 5; 3; 2; 1]; weights = weights/sum(weights);
% weights = ones(7,1)/7;
for ti = 1:npoints
    if ti == 1
        before1 = npoints - 1;
        before2 = npoints - 2;
        before3 = npoints - 3;
    elseif ti == 2
        before1 = 1;
        before2 = npoints - 1;
        before3 = npoints - 2;
    elseif ti == 3
        before1 = 2;
        before2 = 1;
        before3 = npoints - 1;
    else
        before1 = ti - 1;
        before2 = ti - 2;
        before3 = ti - 3;
    end
    
    if ti == npoints
        after1 = 2;
        after2 = 3;
        after3 = 4;
    elseif ti == npoints - 1
        after1 = npoints;
        after2 = 2;
        after3 = 3;
    elseif ti == npoints - 2
        after1 = npoints - 1;
        after2 = npoints;
        after3 = 2;
    else
        after1 = ti + 1;
        after2 = ti + 2;
        after3 = ti + 3;
    end
    
    idxs = [before3, before2, before1, ti, after1, after2, after3];
    vecs(ti,:) = sum(scaled(idxs,:).*weights,1);
    line([b(ti,1); b(ti,1) + vecs(ti,1)],[b(ti,2); b(ti,2) + vecs(ti,2)],'Color','g')
end

axis equal






