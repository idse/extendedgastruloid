function newseg = separate_fused(seg, tau1, tau2, opts)
% take an input segmentation an parameter values and perform approximate
% convex decomposition on each cell in the segmentation; cuts made to
% create approximately convex components are changed from white to black in
% the segmentation

%option to ignore objects under a given area
if ~isfield(opts,'useMinArea')
    opts.useMinArea = false;
end
if opts.useMinArea && ~isfield(opts, 'minArea')
    opts.useMinArea = false;
    disp('No minimum area specified, no minimum area threshold will be applied')
end
if opts.useMinArea
    disp(strcat("Ignoring objects of area less than ", num2str(opts.minArea)))
end
%if troubleshooting option is not set, use default of false
if ~isfield(opts,'flag')
    opts.flag = false;
end
%if diagnostic option is not set, use default of false
if ~isfield(opts,'diagnostics')
    opts.diagnostics = false;
end
%option to ignore holes in objects
if ~isfield(opts,'ignoreholes')
    opts.ignoreholes = false;
end

%if outputing cases closest to the threshold above and below for parameter
%tuning:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%relative concavities of cases nearest to threshold
diagrels = [Inf(1,5); zeros(1,5)];
%absolute concavities of these cases
diagabs = NaN(2,5,2);
%boundaries of relevant objects
diagcoords = cell(2,5);
%cutpoints
cpoints = cell(2,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get boundaries and dependencies from the image
[B, ~, ~, A] = bwboundaries(seg);
n = length(B);
bad = zeros(n,1);

%simplify boundary contours and prevent self-intersections
for bi = 1:n
    b = B{bi};
    if size(unique(b,'rows'),1) > 3 && ~pointsAreCollinear(b)
        B{bi} = correctContours(b);
    else
        bad(bi) = 1;
    end
end
good = bad == 0;
indices = (1:n)';
%holes are objects nested inside of other objects
holes = find(sum(A,2) > 0);
%parents (cells) are boundaries that are not children of other objects
parents = indices(~ismember(indices,holes));
bad = bad(parents);

%organize object contours into a cell array in which the left column has
%arrays of boundary coordinates for parent object and the right column has
%arrays of boundary coordinates for holes (also in a cell at each entry)
cells = cell(length(parents),2);
for ci = 1:length(parents)
    cells{ci,1} = B{parents(ci)};
    if ~opts.ignoreholes
        hidxs = A(:,parents(ci)) > 0;
        hidxs = hidxs & good;
        cells{ci,2} = B(hidxs);
    else
        cells{ci,2} = {};
    end
end
concavecell = cell(size(cells,1),1);
curvecell = concavecell;
veccell = concavecell;

criterion = true;
allcuts = [];
count = 1;
flag = opts.flag;
if flag
    disp(strcat("Number of cells to resolve = ", num2str(size(cells,1))))
end
reps = 0;
while criterion
    reps = reps + 1;
    newcells = cell(1,2);
    newconcave = cell(1);
    newcurve = cell(1);
    newvec = cell(1);
    cuts = [];
    if flag
        disp(strcat("Iteration ", num2str(count)))
    end
    
    for bi = 1:size(cells,1)
        b = cells{bi,1};
        concave = concavecell{bi};
        curve = curvecell{bi};
        vec = veccell{bi};
        
        children = cells{bi,2};
        %get rid of holes with only one or two points
        if ~isempty(children)
            sizes = cellfun(@(x) size(x,1),children);
            children = children(sizes > 5);
        end
        
        %only worry about an object if it has three or more unique vertices
        %and if not all vertices are collinear
        if count == 1
            test = ~bad(bi);
        elseif size(unique(b,'rows'),1) > 3 && ~pointsAreCollinear(b)
            test = true;
        else
            test = false;
        end
        
        if opts.useMinArea && polyarea(b(:,1),b(:,2)) < opts.minArea
            test = false;
        end
        
        if test
            [components, cut, concaves, curves, vecs, diags] = ...
                make_cuts_v2(b, children, concave, curve, vec, tau1, tau2, opts);
            if size(components,1) == 2
                %cut was made
                if opts.diagnostics
                    if diags.relconcave > tau1 && any(diags.relconcave < diagrels(1,:))
                        [~, I] = max(diagrels(1,:));
                        diagrels(1,I) = diags.relconcave;
                        diagabs(1,I,:) = diags.absconcaves;
                        diagcoords{1,I} = b;
                        cpoints{1,I} = diags.cutpoints;
                    end
                end
                newcells = [newcells; components];
                newconcave = [newconcave; concaves];
                newcurve = [newcurve; curves];
                newvec = [newvec; vecs];
            else
                %cut was not made
                if opts.diagnostics
                    if diags.relconcave < tau1 && any(diags.relconcave > diagrels(2,:)) &&...
                        all(diags.absconcaves > tau2)
                        [~, I] = min(diagrels(2,:));
                        diagrels(2,I) = diags.relconcave;
                        diagabs(2,I,:) = diags.absconcaves;
                        diagcoords{2,I} = b;
                        cpoints{2,I} = diags.cutpoints;
                    end
                end                
            end
            cuts = [cuts; cut];
            
        end
        
    end
    if reps > 10000
        error('why is the while loop not terminating')
    end

    allcuts = [allcuts; cuts];
    count = count + 1;
    if size(newcells,1) > 1
        cells = newcells(2:end,:);
        concavecell = newconcave(2:end);
        curvecell = newcurve(2:end);
        veccell = newvec(2:end);
    else
        criterion = false;
    end
end

% apply the cuts by changing pixels along those line segments from white
% to black
imsize = size(seg);
for li = 1:(size(allcuts,1)/2)
    X0 = allcuts(2*li-1,1);
    X1 = allcuts(2*li,1);
    Y0 = allcuts(2*li-1,2);
    Y1 = allcuts(2*li,2);
    if X1 - X0 > Y1 - Y0
        X2 = X0; X3 = X1;
        X4 = X0; X5 = X1;
        Y2 = Y0 + 1; Y3 = Y1 + 1;
        Y4 = Y0 - 1; Y5 = Y1 - 1;
    else
        Y2 = Y0; Y3 = Y1;
        Y4 = Y0; Y5 = Y1;
        X2 = X0 + 1; X3 = X1 + 1;
        X4 = X0 - 1; X5 = X1 - 1;
    end
    xs = [X0, X1; X2, X3; X4, X5];
    xs = max(xs,1); xs = min(xs,imsize(1));
    ys = [Y0, Y1; Y2, Y3; Y4, Y5];
    ys = max(ys,1); ys = min(ys,imsize(2));
    for idx = 1:3
        ns = 0:(1/round(sqrt((xs(idx,2)-xs(idx,1))^2 + (ys(idx,2)-ys(idx,1))^2))):1;
        xn = round(xs(idx,1) +(xs(idx,2) - xs(idx,1))*ns)';
        yn = round(ys(idx,1) +(ys(idx,2) - ys(idx,1))*ns)';
        inds = sub2ind(imsize, xn, yn);
        seg(inds) = 0;
    end
end

%if troubleshooting/optimizing parameters for a given dataset, show closest
%objects above and below the cutoff
if opts.diagnostics
    if isfield(opts,'img')
        img = opts.img;
    else
        img = zeros(size(seg));
    end
    figure
    spc = 1;
    for jidx = 1:2
        for idx = 1:5
            coords = diagcoords{jidx,idx};
            if ~isempty(coords)
                subplot(2,5,spc)
                xl = [min(coords(:,1)),max(coords(:,1))];
                yl = [min(coords(:,2)),max(coords(:,2))];
                height = xl(2) - xl(1);
                width = yl(2) - yl(1);
                xl = xl + [-height/2, height/2];
                yl = yl + [-width/2, width/2];
                imshow(img)
                line(coords(:,1),coords(:,2),'Color','b')
                line(cpoints{jidx,idx}(:,1),cpoints{jidx,idx}(:,2),'Color','r',...
                    'Marker','o');
                title({strcat("Rel concave = ", num2str(diagrels(jidx,idx))),...
                    strcat("Abs concaves = ",...
                    num2str(diagabs(jidx,idx,1)),", ",num2str(diagabs(jidx,idx,2)))})
                if jidx == 1
                    xlabel("Cut")
                else
                    xlabel("Not cut")
                end
                xlim(xl)
                ylim(yl)
            end
            spc = spc + 1;
        end
    end
    set(gcf,'WindowState','maximized')
end


disp(strcat("Number of cuts = ", num2str(size(allcuts,1)/2)))

newseg = seg;


end