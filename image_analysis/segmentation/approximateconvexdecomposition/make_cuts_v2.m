function [components, cuts, concaves, curves, vecss, diagnostics] = ...
    make_cuts_v2(coords, holes, concave, curve, vecs, tau1, tau2, opts)
% perform approximate convex decomposition of an object containing holes
% coords is an ordered n x 2 matrix of x and y values for the polygon
% boundary
% holes is a cell with entries defining the coordinates of boundaries for
% holes inside of the polygon defined by coords
% tau1 is the maximal allowable relative concavity for the approximate 
% decomposition
% tau2 is the maximum allowable absolute concavity for the decomposition
% concave and curve contain precomputed concavities and curvatures from
% previous iterations (if applicable)

%%%add an option not to check for intersections to speed up computation
%%%time

%%%fix the way curvature is computed-> compute it more points out than it
%%%currently uses, and use Lines variable for simpler indexing of points
%%%before and after the selected one; go out as many points as to the
%%%nearest pocket boundary in either direction, up to 5 points

%%%probably would help select the best point to cut from if concavity was
%%%measured as shortest-path instead of straight-line but (a) this will
%%%slow down computation, which is already an issue with the large
%%%micropattern (but less so when carrying over concavities between
%%%iterations) and (b) I'm still not sure how to implement it

%%%use distance to farthest other point in the hole as the absolute
%%%concavity of points in holes? for micropatterns or very dense culture
%%%this becomes extremely computationally expensive but is likely a better
%%%measure (and could be made faster than previous implementation with
%%%pdist function)

%troubleshoot/visualize results at each iteration if flag
flag = opts.flag;
%show over image if one is provided
if isfield(opts,'img')
    img = opts.img;
    flag2 = true;
else
    flag2 = false;
end
%if both cut points exceed some (large) absolute concavity threshold, go
%ahead and cut (potentially a problem for micropattern though?)
if isfield(opts,'tau3')
    tau3 = opts.tau3;
else
    tau3 = 20;
end

npoints = size(coords,1);
%K holds indices of points on the convex hull
K = convhull(coords(:,1), coords(:,2));

%if no existing concavities are provided, all of them should be computed
if isempty(concave)
    concave = NaN(npoints,1);
end
%if no existing curvatures are provided, all of them should be computed
if isempty(curve)
    curve = NaN(npoints,1);
end

if isempty(vecs)
    vecs = NaN(npoints,2);
end

if size(coords,1) ~= size(concave,1)
    error('Check this')
end

%give each pocket on the contour its own index
pockets = zeros(npoints,1);
sorted = sort(unique(K));
pocketidx = 1;
for ki = 1:length(sorted)-1
    if sorted(ki+1) - sorted(ki) > 1
        pockets((sorted(ki)+1):(sorted(ki+1)-1)) = pocketidx;
        pocketidx = pocketidx + 1;
    end
end
%deal with wraparound pocket
if sorted(end) < npoints
    pockets((sorted(end)+1:end)) = pocketidx;
end
if sorted(1) > 1
    pockets(1:sorted(1)-1) = pocketidx;
end

%points on convex hull get (arbitrary) angle measure of pi
curve(K) = pi;
%concavities on convex hull are zero
concave(K) = 0;
%directions on convex hull should probably not be arbitrary but they are
%for now
vecs(K,1) = 1;
vecs(K,2) = 0;

%do cutting: if a cut is made to a hole, incorporate the hull into the
%boundary of the object and calculate corresponding concavities/curvatures
%and if a cut is made to the external boundary, output the newly cut
%boundaries and cut coordinates from the function and stop iterating
%if a potential cut was rejected, stop iterating and output the current
%boundary along with cut coordinates
criterion = true;
reps = 0;
cuts = [];
while criterion
    reps = reps + 1;
%     if reps > 1
%         fprintf('.')
%     end
%     if mod(reps,40) == 0
%         fprintf('\n')
%     end
    
    %get number of holes on each iteration
    nholes = length(holes);
    %make an array of coordinates of all points on boundary or on hole
    %boundaries
    allcoords = cell2mat([coords;holes(:)]);
    
    %lines keeps track of which vertices link to which in polygon formation
    Lines = zeros(size(allcoords));
    Lines(1:npoints,:) = [(1:npoints-1)',(2:npoints)';npoints,1];
    pointidxs = cell(nholes+1,1);
    pointidxs{1} = [zeros(npoints,1), (1:npoints)'];
    indx = npoints + 1;
    for hi = 1:nholes
        nhpoints = size(holes{hi},1);
        pointidxs{hi+1} = [hi*ones(nhpoints,1), (1:nhpoints)'];
        Lines(indx:indx+nhpoints-1,:) = [(indx:indx+nhpoints-2)',...
            (indx+1:indx+nhpoints-1)';indx+nhpoints-1,indx];
        indx = indx + nhpoints;
    end
    pointidxs = cell2mat(pointidxs);
    
    if size(allcoords,1) ~= size(pointidxs,1)
        error('The indexing is messed up somewhere')
    end
    
    %points needing concavity/curvature to be computed are listed in kprime
    kprime = find(isnan(concave));

    %define concavity for each exterior point
    for pj = 1:length(kprime)
        %get boundary points
        pidx = kprime(pj);
        x0 = coords(pidx,1);
        y0 = coords(pidx,2);
        %it may be faster to find pockets first and then iterate over the
        %points in those pockets instead of figuring out the boundaries of the
        %pocket for each individual point

        %k1 is the point on the convex hull preceding the current pocket
        k1 = max(K(K < pidx));
        if isempty(k1)
            k1 = max(K);
        end

        x1 = coords(k1,1);
        y1 = coords(k1,2);
        %k2 is the point on the convex hull following the current pocket
        k2 = min(K(K > pidx));
        if isempty(k2)
            k2 = min(K);
        end
        x2 = coords(k2,1);
        y2 = coords(k2,2);

        %measure straight-line distance to convex hull as concavity
        v = [x2 - x1; y2 - y1];
        u1 = [x0 - x1; y0 - y1];
        u2 = [x0 - x2; y0 - y2];
        cos1 = v'*u1/(norm(v)*norm(u1));
        cos2 = (-v)'*u2/(norm(v)*norm(u2));
        %if the point is not perpendicularly visible to the line segment
        %connecting k1 and k2, use the minimum distance to either k1 or k2 as
        %concavity value
        if cos1 < 0 || cos2 < 0
            concave(pidx) = min(norm(u1), norm(u2));
        else
        %otherwise use perpendicular distance of the point to the line segment
        %between k1 and k2 as concavity
            concave(pidx) = abs(v(2)*x0 - v(1)*y0 + x2*y1 - y2*x1)/norm(v);
        end

        %measure local angle/curvature at each point
        %indexing could be simplified by using Lines to keep track of previous
        %and following points
        %it is probably more robust to use vertices farther out (~5 vertices
        %away) for angle measurement, ensuring that these vertices are still
        %within the pocket
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
        %set curvature/angle measure
        [dot, I] = max([bvec1*avec1',bvec2*avec2']);
        curve(pidx) = acos(dot);
        if orientation(coords(before1,:),coords(pidx,:),coords(after1,:)) == 1
            %if points are oriented clockwise, the angle is towards the
            %convex hull, not inwards towards the polygon interior
            curve(pidx) = 2*pi - curve(pidx);
%             disp('corrected curvature')
        end
        %set normal direction at this point
        if I == 1
            vecs(pidx,:) = -avec1/2 - bvec1/2;
        elseif I == 2
            vecs(pidx,:) = -avec2/2 - bvec2/2;
        end
        
        if norm(vecs(pidx,:)) == 0
            %if the average of the points before and after this one is the
            %same as the coordinates of this point, instead rotate a vector
            %in the direction of the contour by 90 degrees clockwise
            vecs(pidx,:) = ([0,1;-1,0]*(coords(after1,:) - coords(before1,:))')';
        end
        %normalize
        vecs(pidx,:) = vecs(pidx,:)/norm(vecs(pidx,:));
    end
    
    %find the first cut point and its concavity
    crit3 = true;
    sreps = 0;
    mfun = concave./(1+curve);
    while crit3
        sreps = sreps + 1;
        [~, r1] = max(mfun);
        val = concave(r1);
        %remove disallowable points for cutting to
        allowed = ones(npoints,1);
        %don't cut to another point in the same pocket
        pidx = pockets(r1);
        allowed(pockets == pidx) = 0;
        %don't cut to points with absolute concavity of zero
        allowed(concave == 0) = 0;
        %all hole points are allowed
        allowed = [allowed; ones(size(allcoords,1)-npoints,1)]; %ignore warning
        ncandidates = sum(allowed);
        %take the best cut as the shortest one out of the potential candidates, or
        crit2 = true;
        cutreps = 0;
        %pick a cut point, check for intersections with nearby contours, and if
        %there are any try again with a different cut point
        while crit2
            cutreps = cutreps + 1;
            diffcoords = allcoords - coords(r1,:);
            dists = sum(diffcoords.^2,2);
            %normalized vectors in direction of all other points
            dirs = diffcoords./(sqrt(sum(diffcoords.^2,2)));
            vec = vecs(r1,:);
            %dot product of each of these directions with the direction vector
            %of the first selected concave point gives a similarity in
            %direction between 1 (parallel) and -1 (antiparallel)
            dirweight = 1.1 + dirs*vec';%1./(1.1 + dirs*vec');
            %get rid of NaN at r1 because of trying to normalize zero vector
            dirweight(r1) = 1e-3;%Inf;
            %cut lengths to disallowed points are set to infinity
            dists(allowed == 0) = Inf;
            %pick best point for given criteria
            if nholes > 0
                f = (dirweight.^2)./(1+dists);
            elseif nholes == 0
                f = concave./((1 + curve).*(1 + dists)).*dirweight.^2;
            else
                error('Negative holes?')
            end
            f(allowed == 0) = 0;
%             [~, allidx] = min(dists.*dirweight.^2);
            [~,allidx] = max(f);
            cutlen = sqrt(dists(allidx)); %cut length
            hidx = pointidxs(allidx,1); %index of hole (0 if not in a hole)
            r2 = pointidxs(allidx,2); %index of point in hole (r2=allidx if not in a hole)
            %find closest point in original pocket to r2 if r2 is a hole point
%             if hidx > 0
%                 idxs = find(pockets == pidx);
%                 pcoords = coords(idxs,:);
%                 [~,r1] = min(sum((pcoords - allcoords(allidx,:)).^2,2));
%                 r1 = idxs(r1);
%                 val = concave(r1);
%             end

            %check for intersections with lines between nearby vertices
            cut = [coords(r1,:); allcoords(allidx,:)];
            %this search region could be narrowed down to be the intersection
            %of points within the radius of both vertices instead of one or
            %the other
            tests = find(min(sum((allcoords - cut(1,:)).^2,2),...
                sum((allcoords - cut(2,:)).^2,2)) < cutlen.^2);
            excluded = [r1,allidx,Lines(r1,2),Lines(allidx,2),...
                Lines(Lines(:,2) == r1,1)', Lines(Lines(:,2) == allidx,1)'];
            tests = tests(~ismember(tests,excluded));
            intersects = false;
            
            for sidx = 1:length(tests)
                testidx = tests(sidx);
                segment = [allcoords(testidx,:); allcoords(Lines(testidx,2),:)];
                if doIntersect(cut(1,:),cut(2,:),segment(1,:),segment(2,:)) &&...
                        ~any(ismember(cut,segment,'rows'))
                    intersects = true;
                    break
                end
            end

            if intersects
                allowed(allidx) = 0;
            else
                crit2 = false;
                crit3 = false;
            end
            
            if cutreps >= min(100,ncandidates-1)
%                 if reps > 1
%                     fprintf('\n')
%                 end
%                 disp('Trying a new initial cut point')
                mfun(r1) = 0;
                crit2 = false;
            end
            
            if val < tau2 || ncandidates == 0
                %if the first point is not concave enough anyway stop
                %wasting time trying to find a good non-conflicting second
                %cut point
                crit2 = false;
                crit3 = false;
            end
            if ncandidates == 0
                val = 0;
            end
        end
        if sreps > 1000
            error('what now')
        end
    end
    
    if hidx > 0
        b = holes{hidx};
        holecut = true;
        nhpoints = size(b,1);
        %measure curvature at this point in the hole
        %would it be feasible to measure curvature & direction at hole
        %points and also use that for second cut point selection?
        if r2 == nhpoints
            after1 = 2;
            after2 = 3;
        elseif r2 == nhpoints-1
            after1 = nhpoints;
            after2 = 2;
        else
            after1 = r2 + 1;
            after2 = r2 + 2;
        end

        if r2 == 1
            before1 = nhpoints - 1;
            before2 = nhpoints - 2;
        elseif r2 == 2
            before1 = 1;
            before2 = nhpoints - 1;
        else
            before1 = r2 - 1;
            before2 = r2 - 2;
        end
        bvec1 = b(before1,:) - b(r2,:); bvec1 = bvec1/norm(bvec1);
        bvec2 = b(before2,:) - b(r2,:); bvec2 = bvec2/norm(bvec2);
        avec1 = b(after1,:) - b(r2,:); avec1 = avec1/norm(avec1);
        avec2 = b(after2,:) - b(r2,:); avec2 = avec2/norm(avec2);
        holecurve = acos(max(bvec1*avec1',bvec2*avec2'));
    else
        b = coords;
        holecut = false;
    end
    
    if holecut
        testval = 2*val/cutlen;
        %if the cut is in a hole and meets the thresholds cut to the hole, 
        %calculate new concavities for the points added to the object 
        %boundary, and keep looking for potential cuts
        if val > tau2 && (testval > tau1 || (sum([curve(r1),holecurve] < pi/2) > 0))
            nhpoints = size(b,1);
            %determine two new points, half a pixel apart from the initial
            %point for each of these            
            v1 = coords(Lines(:,2) == r1,:) - coords(r1,:);
            v1 = 0.5*v1/max(1,norm(v1));
            v2 = coords(Lines(r1,2),:) - coords(r1,:);
            v2 = 0.5*v2/max(1,norm(v2));
            u1 = allcoords(Lines(:,2) == allidx,:) - allcoords(allidx,:);
            u1 = 0.5*u1/max(1,norm(u1));
            u2 = allcoords(Lines(allidx,2),:) - allcoords(allidx,:);
            u2 = 0.5*u2/max(1,norm(u2));

            coord1 = [r1+1:npoints, 1:r1-1];
            b1 = [r2-1:-1:1, nhpoints:-1:r2+1];
            components = cell(1,2);
            
            newcoords = [...
                coords(r1,:) + v2;...
                coords(coord1,:);...
                coords(r1,:) + v1;...
                b(r2,:) + u1;...
                b(b1,:);...
                b(r2,:) + u2];
            
            cut = [coords(r1,:); b(r2,:)];
            cuts = [cuts; cut];
            idxs = 1:length(holes);
            idxs = idxs(idxs ~= hidx);
            newholes = holes(idxs);
        else
            criterion = false;
            %if only one component is returned by the function, that
            %indicates that cutting this contour is done
            components = cell(1,2);
            components{1,1} = coords;
            components{1,2} = holes;
            concaves = [];
            curves = [];
            vecss = [];
            
            diagnostics = struct(...
                'cutpoints',    [coords(r1,:); b(r2,:)],...
                'relconcave',   testval,...
                'absconcaves',  [val, 0]);
        end
        
    else
        %if the cut is to another boundary point and meets the threshold,
        %cut the object into two pieces, save the cut and exit the function
        testval = (val + concave(r2))/cutlen;
        if min(val,concave(r2)) > tau2 && (testval > tau1 || (sum(curve([r1,r2]) < pi/2) == 2) || (concave(r1) > tau3 && concave(r2) > tau3))            
            components = cell(2,2);
            concaves = cell(2,1);
            curves = cell(2,1);
            vecss = cell(2,1);
            newconcave = concave;
            newcurve = curve;
            newvecs = vecs;
            %make curvature and concavity values inside the two pockets
            %that are cut between be NaN to indicate that they should be
            %recomputed at the next iteration
            nanidxs = ismember(pockets,[pockets(r1),pockets(r2)]);
            newconcave(nanidxs) = NaN;
            newcurve(nanidxs) = NaN;
            newvecs(nanidxs,:) = NaN;
            
            %coord1 and coord2 are indices of the boundary vertices that go
            %into each of the new contours, listed in order and with the
            %first point repeated as the last so that contours are closed
            coord1 = min(r1,r2):max(r1,r2);
            coord2 = [1:min(r1,r2),max(r1,r2):npoints];
            if numel(coord1) < 3 || numel(coord2) < 3
                error('bad cut')
            end
            components{1,1} = coords(coord1,:);
            components{2,1} = coords(coord2,:);
            concaves{1} = newconcave(coord1);
            concaves{2} = newconcave(coord2);
            curves{1} = newcurve(coord1);
            curves{2} = newcurve(coord2);
            vecss{1} = newvecs(coord1,:);
            vecss{2} = newvecs(coord2,:);
            cuts = [cuts; cut];
            
            assignments = zeros(length(holes),1);
            for hi = 1:length(holes)
                point = holes{hi}(1,:);
                if inpolygon(point(1),point(2),components{1,1}(:,1),components{1,1}(:,2))
                    assignments(hi) = 1;
                else
                    assignments(hi) = 2;
                end
            end
            components{1,2} = holes(assignments == 1);
            components{2,2} = holes(assignments == 2);
        else
            %if the cut is to another boundary point and does not meet the
            %criteria, output the original boundary with holes included
            concaves = {[]};
            curves = {[]};
            vecss = {[]};
            components = cell(1,2);
            components{1,1} = coords;
            components{1,2} = holes;
        end
        criterion = false;
        diagnostics = struct(...
            'cutpoints',    coords([r1,r2],:),...
            'relconcave',   testval,...
            'absconcaves',  [val, concave(r2)]);
    end
    
    %visualize cutting results if the object is cut in half (cut
    %successfully made from exterior boundary to exterior boundary)
    if flag && criterion == false% ~holecut && (size(components,1) == 2)
        xl = [min(coords(:,1)),max(coords(:,1))];
        yl = [min(coords(:,2)),max(coords(:,2))];
        height = xl(2) - xl(1);
        width = yl(2) - yl(1);
        xl = xl + [-height/2, height/2];
        yl = yl + [-width/2, width/2];
        
        cla
        if flag2
            imshow(img)
            hold on
        end
        
        if size(cuts,1) > 1
            cut1 = reshape(cuts(:,1),2,[]);
            cut2 = reshape(cuts(:,2),2,[]);
            line(cut1,cut2,'Marker','.','Color','r')
        end
        comp1 = components{1,1};
        scatter(allcoords([r1,allidx],1),allcoords([r1,allidx],2),...
            36,'m','filled')
        line(coords(K,1), coords(K,2), 'Color', 'g')
        line(comp1(:,1), comp1(:,2),'Color','c')
        if size(components,1) == 2
            comp2 = components{2,1};
            line(comp2(:,1), comp2(:,2),'Color','b')
            for ti = 1:length(components{2,2})
                hole = components{2,2}{ti};
                line(hole(:,1), hole(:,2),'Color','m')
            end
            title(strcat("Cut into parts; Abs concaves = ", num2str(val),...
                ", ", num2str(concave(r2)), "Rel concave = ",...
                num2str(testval)))
        else
            title(strcat("Not cut; Abs concaves = ", num2str(val),", ",...
                num2str(concave(r2)),", Rel concave = ", num2str(testval)))
        end
        
        for ti = 1:length(components{1,2})
            hole = components{1,2}{ti};
            line(hole(:,1), hole(:,2),'Color','m')
        end
        xlim(xl)
        ylim(yl)
        xlabel(strcat("Curvature 1 = ", num2str(curve(r1)*180/pi),...
            "; Curvature2 = ", num2str(curve(r2)*180/pi)))
        pause
    end
    
    if criterion && holecut
        %move previously computed features into new sets of indices for
        %expanded boundary
        newconcave = [NaN;...
            concave(coord1);...
            NaN(length(b1)+3,1)];
        newcurve = [NaN;...
            curve(coord1);...
            NaN(length(b1)+3,1)];
        newvecs = [NaN(1,2);...
            vecs(coord1,:);...
            NaN(length(b1)+3,2)];
%         cidxs = find(isnan(newconcave));
        thing = zeros(npoints,1);
        thing(K) = 1;
        thing = [NaN;...
            thing(coord1);...
            NaN(length(b1)+3,1)];
        K = find(thing == 1);
        pocketidx = pockets(r1);
        newpockets = [pocketidx;...
            pockets(coord1);...
            pocketidx*ones(length(b1)+3,1)];
        holes = newholes;
        
        concave = newconcave;
        curve = newcurve;
        vecs = newvecs;
        pockets = newpockets;
        coords = newcoords;
        npoints = size(coords,1);        
    end
    

end
% if reps > 1
%     fprintf('\n')
%     disp(strcat("Number of iterations = ", num2str(reps)))
% end


end