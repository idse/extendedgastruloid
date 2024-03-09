function newcoords = correctContours(coords)
%function to resolve self-intersecting and overlapping object contours by
%finding and separating repeated points
%repeats are handled in the order in which they occur, and resolved such
%that coordinates are 

%this function is well-behaved for contours output from bwboundaries, but
%not necessarily for generalized or modified contours
coords = coords(1:end-1,:);
npoints = size(coords,1);
Lines = [(1:npoints)',[2:npoints,1]'];
[~, ia, ic] = unique(coords,'rows','stable');
%repeated points are the ones in the original array but not the unique one
repeats = setdiff((1:npoints)',ia);
newcoords = coords;
flag = false;
badidx = zeros(size(newcoords,1),1);

for ri = 1:size(repeats,1)
    %indices of repeated points
    idxs = find(ic == ic(repeats(ri)));
    %coordinates of the repeated point
    xy = coords(idxs(1),:);
    nidxs = size(idxs,1);
    %after points
    af = zeros(nidxs,2);
    %before points
    bf = zeros(nidxs,2);
    us = zeros(nidxs,2);
    thetas = zeros(nidxs,1);
    
    for ti = 1:nidxs
        bf(ti,:) = coords(Lines(:,2) == idxs(ti),:);
        af(ti,:) = coords(Lines(idxs(ti),2),:);
    end
    
    for ti = 1:nidxs
        if ti == 1
            us(ti,:) = mean([bf(ti,:);af(nidxs,:)],1) - xy;
            if norm(us(ti,:)) == 0
                us(ti,:) = bf(ti,:) - xy;
            end
        else
            us(ti,:) = mean([bf(ti,:); af(ti-1,:)],1) - xy;
            if norm(us(ti,:)) == 0
                us(ti,:) = bf(ti,:) - xy;
            end
        end
        us(ti,:) = 0.5*us(ti,:)/norm(us(ti,:));
    end
    
    for ti = 1:nidxs
        if ti == nidxs
            thetas(ti) = acos(us(1,:)*us(nidxs,:)')/2;
        else
            thetas(ti) = acos(us(ti,:)*us(ti+1,:)')/2;
        end
        newcoords(idxs(ti),:) = xy +...
            ([cos(thetas(ti)) sin(thetas(ti)); -sin(thetas(ti)) cos(thetas(ti))]*us(ti,:)')';
        if sum(isnan(newcoords(idxs(ti),:)),'all') > 0
            disp('why are there nans')
            badidx(idxs(ti)) = 1;
            flag = true;
        end
    end
end

if flag
    newcoords(badidx == 1,1) = xy(1);
    newcoords(badidx == 1,2) = xy(2);
    clf
    line(newcoords(:,1),newcoords(:,2))
    pause
end

% newcoords = [newcoords; newcoords(1,:)];


end