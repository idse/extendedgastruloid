function [components, cut, diagnostics] = ...
    approx_decomp_v2(coords, tau1, tau2, s, t, opts)
% coords is an ordered n x 2 matrix of x and y values for the polygon boundary
% tau is the maximal allowable concavity for the approximate decomposition
% s and t are parameters that tune the way in which cuts are selected;
% higher s selects for minimizing concavity, and higher t selects for
% shorter cuts
npoints = size(coords,1);
K = convhull(coords(:,1), coords(:,2));
kprime = find(~ismember(1:npoints,K))';
concave = zeros(npoints,1);
curves = 180*ones(npoints,1); %points on 
flag = opts.flag;
if isfield(opts,'img')
    img = opts.img;
    flag2 = true;
else
    flag2 = false;
end
if ~isfield(opts,'flag3')
    flag3 = false;
else
    flag3 = opts.flag3;
end

area = polyarea(coords(:,1),coords(:,2));

for pj = 1:length(kprime)
    %get boundary points
    pidx = kprime(pj);
    x0 = coords(pidx,1);
    y0 = coords(pidx,2);
    
    k1 = max(K(K < pidx));
    if isempty(k1)
        k1 = max(K);
    end
    x1 = coords(k1,1);
    y1 = coords(k1,2);
    
    k2 = min(K(K > pidx));
    if isempty(k2)
        k2 = min(K);
    end
    x2 = coords(k2,1);
    y2 = coords(k2,2);
    
    v = [x2 - x1; y2 - y1];
    u1 = [x0 - x1; y0 - y1];
    u2 = [x0 - x2; y0 - y2];
    cos1 = v'*u1/(norm(v)*norm(u1));
    cos2 = (-v)'*u2/(norm(v)*norm(u2));
    
    if cos1 < 0 || cos2 < 0
        concave(pidx) = min(norm(u1), norm(u2));
    else
        concave(pidx) = abs(v(2)*x0 - v(1)*y0 + x2*y1 - y2*x1)/norm(v);
    end
    
    %measure local angle/curvature at each point
    if pidx == npoints
        after1 = 2;
        after2 = 3;
    elseif pidx == npoints-1
        after1 = npoints;
        after2 = 2;
    else
        after1 = pidx + 1;
        after2 = pidx + 2;
    end

    if pidx == 1
        before1 = npoints - 1;
        before2 = npoints - 2;
    elseif pidx == 2
        before1 = 1;
        before2 = npoints - 1;
    else
        before1 = pidx - 1;
        before2 = pidx - 2;
    end
    bvec1 = coords(before1,:) - coords(pidx,:); bvec1 = bvec1/norm(bvec1);
    bvec2 = coords(before2,:) - coords(pidx,:); bvec2 = bvec2/norm(bvec2);
    avec1 = coords(after1,:) - coords(pidx,:); avec1 = avec1/norm(avec1);
    avec2 = coords(after2,:) - coords(pidx,:); avec2 = avec2/norm(avec2);
    curves(pidx) = acos(max(bvec1*avec1',bvec2*avec2'))*180/pi;
end

criterion2 = true;
choosereps = 0;
while criterion2
    choosereps = choosereps + 1;
    if choosereps > npoints - 1
        error('no good point to choose?')
    end
    [~,r] = max(concave./(0.5+curves/180));
    val = concave(r);
    components = cell(2,1);

    k1 = max(K(K < r));
    if isempty(k1)
        k1 = max(K);
    end

    k2 = min(K(K > r));
    if isempty(k2)
        k2 = min(K);
    end
    if k1<k2
        pocket = k1:k2;
    else
        pocket = [k1:npoints,1:k2];
    end

    newcoord = coords;
    newcoord(r,:) = [Inf, Inf];
    f = 1 + s*(concave./curves*180)./(t*sum((newcoord - coords(r,:)).^2,2));
    f(pocket) = 0;
    % need to eliminate points that would not resolve the notch??

    % assign a cut, then check if it intersects the polygon; if it does,
    % mark that link as bad and try again
    criterion = true;
    reps = 0;
    while criterion
        reps = reps + 1;
        [~, x] = max(f);
        count = 0;
        if x == 1 || r == 1
            nocheck = [x,r,x-1,r-1,npoints,npoints-1];
        else
            nocheck = [x,r,x-1,r-1];
        end
        checkpoints = 1:npoints;
        checkpoints = checkpoints(~ismember(checkpoints,nocheck));
        for pnt = 1:length(checkpoints)
            point = checkpoints(pnt);
            if point == npoints
                point2 = 1;
            else
                point2 = point + 1;
            end
            if doIntersect(coords(r,:),coords(x,:),coords(point,:),coords(point2,:))
                count = count + 1;
            end
        end

        if count > 0
            %if the cut intersected at least one side of the polygon,
            %discard this point and try again
            f(x) = 0;
        else
            criterion = false;
            criterion2 = false;
        end
        if reps > npoints
            disp('No good cut point, choosing another initial point')
            criterion = false;
            concave(r) = 0;
        end
    end
end

%compare the concavities of x and r to the cut length
%high concavities and small cut length result
testval = (concave(x) + concave(r))/(norm(coords(x,:) - coords(r,:)));
curves = curves([r,x]);

diagnostics = struct;
diagnostics.concave1 = concave(r);
diagnostics.concave2 = testval;
diagnostics.curve1 = curves(1);
diagnostics.curve2 = curves(2);
diagnostics.area = area;

xl = [min(coords(:,1)),max(coords(:,1))];
yl = [min(coords(:,2)),max(coords(:,2))];
diagnostics.boundingbox = [xl; yl];

height = xl(2) - xl(1);
width = yl(2) - yl(1);
xl = xl + [-height/2, height/2];
yl = yl + [-width/2, width/2];

if (testval < tau1 && (sum(curves < 80) < 2)) || (val < tau2)
    diagnostics.split = false;
    components = {coords};
    cut = [];
    if flag
        cla
        if flag2
            imshow(img)
            hold on
        end
        scatter(coords(r,1), coords(r,2), 36, 'r', 'filled')
        hold on
        scatter(coords(x,1), coords(x,2), 36, 'g', 'filled')
        line(coords(K,1), coords(K,2), 'Color', 'g')
        line(coords(:,1), coords(:,2), 'Color', 'b')
        xlim(xl)
        ylim(yl)
        title(strcat("Relative Concavity = ", num2str(testval),...
            ", Absolute Concavity = ", num2str(val),", Area = ", num2str(area)))
        xlabel(strcat("Curvature 1 = ", num2str(curves(1)),...
            "; Curvature 2 = ", num2str(curves(2))))
        if ~flag3
            pause
        end
    end
else
    diagnostics.split = true;
    components{1} = [coords(min(x,r):max(x,r),:);coords(min(x,r),:)];
    components{2} = [coords(1:min(x,r),:);coords(max(x,r):end,:)];
    cut = [coords(r,:); coords(x,:)];
    
    if flag
        cla
        if flag2
            imshow(img)
            hold on
        end
        line(coords(K,1), coords(K,2), 'Color', 'g')
%         line([coords(r,1);coords(x,1)],[coords(r,2);coords(x,2)],'Color','r')
        line(components{1}(:,1), components{1}(:,2),'Color','c')
        line(components{2}(:,1), components{2}(:,2),'Color','m')
        xlim(xl)
        ylim(yl)
        xlabel(strcat("Curvature 1 = ", num2str(curves(1)),...
            "; Curvature 2 = ", num2str(curves(2)),", Area = ",...
            num2str(area)))
        title(strcat("Relative Concavity = ",...
            num2str(testval), ", Absolute Concavity = ", num2str(val)))
        if ~flag3
            pause
        end
    end
end

if flag && flag3
    f = gcf;
    f.Position(1) = 0;
    criterion = true;
    while criterion
        [indx1, tf] = listdlg('PromptString',{'Should this be split?'},...
            'ListString',...
            {'No','Yes','Bad nucleus','Dividing nucleus','Pause'});

        if indx1 == 5
            pause
        else
            criterion = false;
        end
    end

    if indx1 == 2
        div = false;
        split = true;
        bad = false;
    elseif indx1 == 1
        div = false;
        split = false;
        bad = false;
    elseif indx1 == 3
        div = [];
        split = [];
        bad = true;
    else
        div = true;
        bad = false;
        split = false;
    end

    diagnostics.dividing = div;
    diagnostics.manual = split;
    diagnostics.bad = bad;

    if tf == 0
        error('Force exit')
    end
end

end