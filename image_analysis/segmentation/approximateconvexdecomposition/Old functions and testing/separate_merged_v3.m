function [newseg, diagnostics] =...
    separate_merged_v3(seg, tau1, tau2, s, t, opts)
% take an input segmentation an parameter values and perform approximate
% convex decomposition on each cell in the segmentation; cuts made to
% create approximately convex components are changed from white to black in
% the segmentation

B = bwboundaries(seg);
n = length(B);
bad = zeros(n,1);
for bi = 1:n
    b = B{bi};
    if size(unique(b,'rows'),1) > 3 && ~pointsAreCollinear(b)
        B{bi} = correctContours(b);
    else
        bad(bi) = 1;
    end
end

criterion = true;
allcuts = [];
count = 1;
flag = opts.flag;
if flag
    disp(strcat("Number of cells to resolve = ", num2str(length(B))))
end
if isfield(opts, 'diagnostics')
    diagnostics = opts.diagnostics;
else
    diagnostics = regionprops(seg,'MajorAxisLength','MinorAxisLength');
end

while criterion
    bs = cell(size(B,1),2);
    cuts = cell(size(B));
    if flag
        disp(strcat("Iteration ", num2str(count)))
    end
    
    if count > 1
        opts.flag = false;
    end
    
    for bi = 1:size(B,1)
%         fprintf('.')
        disp(strcat(num2str(bi), " of ", num2str(length(B))))
        b = B{bi};
%         b = b(1:end-1);
        if bad(bi) == 0
            [components, cut, temp] = approx_decomp_v2(b, tau1, tau2, s, t, opts);
        else
            components = {bi};
            cut = [];
        end
        cuts{bi} = cut;
        bs{bi,1} = components{1};
        if numel(components) == 2
            bs{bi,2} = components{2};
        end
        
        if count == 1
            fields = fieldnames(temp);
            for fi = 1:length(fields)
                diagnostics(bi).(fields{fi}) = temp.(fields{fi});
            end
        end
    end
%     fprintf('/n')
    
    %find nuclei that were cut; try decomposing them again
    goods = cellfun(@(x) ~isempty(x), bs(:,2));
    if sum(goods) == 0
        criterion = false;
    end
    B = bs(goods,:);
    B = B(:);
    cuts = cell2mat(cuts(goods));
    allcuts = [allcuts; cuts];
    count = count + 1;
    if count > 1000
        error('why')
    end
end

% disp(strcat("Number of cuts = ",num2str(size(allcuts,1)/2)))

% apply the cuts by changing pixels along those line segments from white
% to black
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
    ys = [Y0, Y1; Y2, Y3; Y4, Y5];
    for idx = 1:3
        ns = 0:(1/round(sqrt((xs(idx,2)-xs(idx,1))^2 + (ys(idx,2)-ys(idx,1))^2))):1;
        xn = round(xs(idx,1) +(xs(idx,2) - xs(idx,1))*ns)';
        yn = round(ys(idx,1) +(ys(idx,2) - ys(idx,1))*ns)';
        inds = sub2ind(size(seg), xn, yn);
        seg(inds) = 0;
    end
end


disp(strcat("Number of cuts = ", num2str(size(allcuts,1)/2)))

newseg = seg;


end