%% Tracking parameters
colony_center = [604.5, 591];
pixel_size = 0.6029203;
last_frame = 146;
first_frame = 1;
tracking_radius = 175 / pixel_size;
%% Loading MIP image and segmentation mask
img = tiffreadVolume('230611_M306_TN_10MIN_B50_HEP_BY_MK_4_37H_FusionStitcher_MIP_p0001_w0000_MIP.tif');
mask = tiffreadVolume('230611_M306_TN_10MIN_B50_HEP_BY_MK_4_37H_FusionStitcher_MIP_p0001_w0000_MIP_cp_masks.tif');
%% Loading tracking result
load('tracking_inward_outward.mat');
mask_out_only = imread('LblImg_mask_out_cell_only.tif');
%% Visualize tracks
img_tracking_frame = img(:, :, last_frame);
img_first_frame = img(:, :, 1);
img_first_frame = imadjust(img_first_frame, [0.06, 0.65], [0, 1], 0.5);
img_tracking_frame = imadjust(img_tracking_frame, [0.1, 0.5], [0, 0.6], 0.5);
C = imfuse(img_tracking_frame, img_first_frame, 'Scaling', 'none');
f = figure;
imshow(C)
bw_circle = createMask(images.roi.Circle(gca, 'Center', colony_center, 'Radius', 570));
close(f);
C(cat(3, ~bw_circle, ~bw_circle, ~bw_circle)) = 255;
size_incre = 0.1;
track_overlay = img_first_frame + img_tracking_frame;
colormap_track = cmap_rand(1e5);
parent_id_plotted = [];
transparency = 0.4;
f = figure;
imshow(C);
hold on
for i = tracked_cell_ids
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    daughter_last_pos = positions(end, :);
    pixel_pos = round(daughter_last_pos);
    for k = 1:height(frames)
        marker_size = (frames(k) - frames(1) + 1) * size_incre;
        plot(positions(k, 1), positions(k, 2), 'Marker', '.', 'Color', [colormap_track(i, :), transparency], 'MarkerSize', min(marker_size, 4));
        hold on
    end
    plot(positions(:, 1), positions(:, 2), 'Color', [colormap_track(i, :), transparency], 'LineWidth', 1);
    hold on
    if height(positions) > 1
        positions = cell_position_tracking(i).positions;
        frames = cell_position_tracking(i).frames;
        daughter_first_pos = positions(1, :);
        parent = cell_position_tracking(i).parent;
        while parent ~= 0 && ~ismember(parent, parent_id_plotted)
            positions = cell_position_tracking(parent).positions;
            frames = cell_position_tracking(parent).frames;
            plot(positions(end, 1), positions(end, 2), 'Marker', 'pentagram', 'Color', [colormap_track(parent, :), transparency], 'MarkerSize', 4, 'MarkerFaceColor', 'none', 'LineWidth', 1.2);
            hold on
            positions = cat(1, positions, daughter_first_pos);
            for k = 1:height(frames)
                marker_size = (frames(k) - frames(1) + 1) * size_incre;
                plot(positions(k, 1), positions(k, 2), 'Marker', '.', 'Color', [colormap_track(parent, :), transparency], 'MarkerSize', min(marker_size, 4));
                hold on
            end
            if height(positions) > 1
                plot(positions(:, 1), positions(:, 2), 'Color', [colormap_track(parent, :), transparency], 'LineWidth', 1);
                hold on
            end
            parent_id_plotted(end + 1, 1) = parent;
            parent = cell_position_tracking(parent).parent;
            daughter_first_pos = positions(1, :);
        end
    end
end
text(0.88, 0.96, sprintf("%.3gh", 36), 'Color', [1, 0, 1], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
text(0.88, 0.92, sprintf("%.3gh", 61), 'Color', [0, 1, 0], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
hold on
ax = gca;
pixPerUm = 1.6586;
scalebarLength = 50; 
unit = sprintf('%sm', '\mu');
hScalebar = scalebar(ax, 'x', scalebarLength, unit, 'Location', 'southwest', 'ConversionFactor', pixPerUm, 'FontSize', 20, 'FontName', 'Arial', 'Margin', 40, 'LineWidth', 2);
exportgraphics(ax, "full_tracks.png", 'Resolution', 300);
close(f)
%% Total displacement vectors for both inward moving cells only
img_tracking_frame = img(:, :, last_frame);
img_first_frame = img(:, :, 1);
img_first_frame = imadjust(img_first_frame, [0.06, 0.65], [0, 1], 0.5);
img_tracking_frame = imadjust(img_tracking_frame, [0.1, 0.5], [0, 0.6], 0.5);
C = imfuse(img_tracking_frame, img_first_frame, 'Scaling', 'none');
f = figure;
imshow(C)
bw_circle = createMask(images.roi.Circle(gca, 'Center', colony_center, 'Radius', 570));
close(f);
C(cat(3, ~bw_circle, ~bw_circle, ~bw_circle)) = 255;
size_incre = 0.1;
track_overlay = img_first_frame + img_tracking_frame;
colormap_track = cmap_rand(1e5);
parent_id_plotted = [];
f = figure;
imshow(C);
hold on
for i = tracked_cell_ids
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    daughter_last_pos = positions(end, :);
    pixel_pos = round(daughter_last_pos);
    hold on
    parent = cell_position_tracking(i).parent;
    while parent ~= 0
        positions = cat(1, cell_position_tracking(parent).positions, positions);
        parent = cell_position_tracking(parent).parent;
    end
    vector_length = norm(daughter_last_pos - positions(1, :));
    if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
        quiver(positions(1, 1), positions(1, 2), daughter_last_pos(1) - positions(1, 1), daughter_last_pos(2) - positions(1, 2), 0, 'Color', 'white', 'LineWidth', 1.2, 'MaxHeadSize', 50/vector_length);
    end
end
text(0.88, 0.96, sprintf("%.3gh", 36), 'Color', [1, 0, 1], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
text(0.88, 0.92, sprintf("%.3gh", 61), 'Color', [0, 1, 0], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
hold on
ax = gca;
pixPerUm = 1.6586;
scalebarLength = 50; 
unit = sprintf('%sm', '\mu'); 
hScalebar = scalebar(ax, 'x', scalebarLength, unit, 'Location', 'southwest', 'ConversionFactor', pixPerUm, 'FontSize', 20, 'FontName', 'Arial', 'Margin', 40, 'LineWidth', 2);
exportgraphics(ax, "mapping_total_displacement_inward_only.png", 'Resolution', 300);
close(f)
%% Total displacement vectors for both inward and outward moving cells
img_tracking_frame = img(:, :, last_frame);
img_first_frame = img(:, :, 1);
img_first_frame = imadjust(img_first_frame, [0.06, 0.65], [0, 1], 0.5);
img_tracking_frame = imadjust(img_tracking_frame, [0.1, 0.5], [0, 0.6], 0.5);
C = imfuse(img_tracking_frame, img_first_frame, 'Scaling', 'none');
f = figure;
imshow(C)
bw_circle = createMask(images.roi.Circle(gca, 'Center', colony_center, 'Radius', 570));
close(f);
C(cat(3, ~bw_circle, ~bw_circle, ~bw_circle)) = 255;
size_incre = 0.1;
track_overlay = img_first_frame + img_tracking_frame;
colormap_track = cmap_rand(1e5);
parent_id_plotted = [];
f = figure;
imshow(C);
hold on
for i = tracked_cell_ids
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    daughter_last_pos = positions(end, :);
    pixel_pos = round(daughter_last_pos);
    hold on
    parent = cell_position_tracking(i).parent;
    while parent ~= 0
        positions = cat(1, cell_position_tracking(parent).positions, positions);
        parent = cell_position_tracking(parent).parent;
    end
    vector_length = norm(daughter_last_pos - positions(1, :));
    if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
        quiver(positions(1, 1), positions(1, 2), daughter_last_pos(1) - positions(1, 1), daughter_last_pos(2) - positions(1, 2), 0, 'Color', 'white', 'LineWidth', 1.2, 'MaxHeadSize', 50/vector_length);
    else
        quiver(positions(1, 1), positions(1, 2), daughter_last_pos(1) - positions(1, 1), daughter_last_pos(2) - positions(1, 2), 0, 'Color', [0, 1, 1], 'LineWidth', 1.2, 'MaxHeadSize', 50/vector_length);
    end
end
text(0.88, 0.96, sprintf("%.3gh", 36), 'Color', [1, 0, 1], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
text(0.88, 0.92, sprintf("%.3gh", 61), 'Color', [0, 1, 0], 'FontSize', 20, 'FontName', 'Arial', 'Units', 'normalized', 'FontWeight','bold');
hold on
ax = gca;
pixPerUm = 1.6586;
scalebarLength = 50; 
unit = sprintf('%sm', '\mu'); 
hScalebar = scalebar(ax, 'x', scalebarLength, unit, 'Location', 'southwest', 'ConversionFactor', pixPerUm, 'FontSize', 20, 'FontName', 'Arial', 'Margin', 40, 'LineWidth', 2);
exportgraphics(ax, "mapping_total_displacement_in_out.png", 'Resolution', 300);
close(f)
%% Radial mapping for inward moving cells only
first_frame_mask = mask(:, :, first_frame);
rp = regionprops("table", first_frame_mask, "Centroid");
rp = rp(~ismissing(rp), :);
radial_dist_first_full = sqrt(sum((rp.Centroid - colony_center).^2, 2)) * pixel_size;
radial_dist_last_full = [];
radial_dist_last = [];
mapped_radial_position = [];
for i = cell_ids_last_frame
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    parent = cell_position_tracking(i).parent;
    daughter_last_pos = positions(end, :);
    radial_dist_last_full(end + 1, 1) = norm(daughter_last_pos - colony_center) * pixel_size;
    if ~ismember(i, tracked_cell_ids)
        continue;
    end
    pixel_pos = round(daughter_last_pos);
    parent = cell_position_tracking(i).parent;
    while parent ~= 0
        positions = cat(1, cell_position_tracking(parent).positions, positions);
        frames = cat(1, cell_position_tracking(parent).frames, frames);
        parent = cell_position_tracking(parent).parent;
    end  
    positions = positions(frames >= first_frame, :);
    frame = frames(frames >= first_frame);
    position_earlist = positions(1, :);
    frame_earliest = frames(1, :);
    if frame_earliest == first_frame
        if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
            mapped_radial_position(end + 1, 1) = norm(position_earlist - colony_center) * pixel_size;
            radial_dist_last(end + 1, 1) = norm(daughter_last_pos - colony_center) * pixel_size;
        end
    end
end


f = figure;
f.Position = [10, 10, 2000, 2000];
[fp,xfp] = kde(mapped_radial_position, 'Support', 'positive');
h(4) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [1 0 1], "FaceColor", [1 0 1]);
hold on
[fp,xfp] = kde(radial_dist_last_full, 'Support', 'positive');
h(3) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [0 0.5 0], "FaceColor", [0 0.5 0]);
[fp,xfp] = kde(radial_dist_last', 'Support', 'positive');
%histogram(radial_dist_last, 'Normalization', 'pdf', 'BinMethod', 'fd', 'FaceColor', [0 0.4470 0.7410]);
h(2) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [0 1 0], "FaceColor", [0 1 0]);
[fp,xfp] = kde(radial_dist_first_full', 'Support', 'positive');
h(1) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [0.5 0 0.5], "FaceColor", [0.5 0 0.5]);
ax = gca;
set(ax, 'FontSize', 34.5);
le = legend([h(1) h(3) h(4) h(2)], ...
    {sprintf('36h (n=%d)', length(radial_dist_first_full)), ...
    sprintf('61h (n=%d)', length(radial_dist_last_full)), ...
    sprintf('center 36h (n=%d)', length(mapped_radial_position)), ...
    sprintf('center 61h (n=%d)', length(radial_dist_last))}, ...
    'Location','northwest', 'FontSize', 26, 'Box','off');
set(ax, "FontWeight", 'bold');
fontname(ax, "Arial");
axis square
set(ax, 'Units', 'centimeters', 'OuterPosition', [0 0 28 28]);
set(ax, "LineWidth", 4)
xlabel(ax, "r (µm)")
yticks([0, 0.02])
ylab = ylabel(ax, "density");
ylab.Position(1) = ylab.Position(1) + 50;
xticks([0, 100, 200, 300, 400]);
%}
box off
exportgraphics(ax,'mapping_radial_distribution_inward_only.png','Resolution',300)
close(f)
%% Radial mapping for both inward and outward moving cells
first_frame_mask = mask(:, :, first_frame);
rp = regionprops("table", first_frame_mask, "Centroid");
rp = rp(~ismissing(rp), :);
radial_dist_first_full = sqrt(sum((rp.Centroid - colony_center).^2, 2)) * pixel_size;
radial_dist_last_full = [];
radial_dist_last_in = [];
radial_dist_last_out = [];
mapped_radial_position_in = [];
mapped_radial_position_out = [];
for i = cell_ids_last_frame
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    parent = cell_position_tracking(i).parent;
    daughter_last_pos = positions(end, :);
    radial_dist_last_full(end + 1, 1) = norm(daughter_last_pos - colony_center) * pixel_size;
    if ~ismember(i, tracked_cell_ids)
        continue;
    end
    pixel_pos = round(daughter_last_pos);
    parent = cell_position_tracking(i).parent;
    while parent ~= 0
        positions = cat(1, cell_position_tracking(parent).positions, positions);
        frames = cat(1, cell_position_tracking(parent).frames, frames);
        parent = cell_position_tracking(parent).parent;
    end  
    positions = positions(frames >= first_frame, :);
    frame = frames(frames >= first_frame);
    position_earlist = positions(1, :);
    frame_earliest = frames(1, :);
    if frame_earliest == first_frame
        if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
            mapped_radial_position_in(end + 1, 1) = norm(position_earlist - colony_center) * pixel_size;
            radial_dist_last_in(end + 1, 1) = norm(daughter_last_pos - colony_center) * pixel_size;
        elseif mask_out_only(pixel_pos(2), pixel_pos(1)) ~= 0
            mapped_radial_position_out(end + 1, 1) = norm(position_earlist - colony_center) * pixel_size;
            radial_dist_last_out(end + 1, 1) = norm(daughter_last_pos - colony_center) * pixel_size;
        end
    end
end
h = [];
f = figure;
f.Position = [10, 10, 2000, 2000];
[fp,xfp] = kde(mapped_radial_position_in, 'Support', 'positive');
h(6) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, 'LineStyle', ':', "EdgeColor", [0.8 0 0.8], "FaceColor", [0.8 0 0.8]);
hold on
[fp,xfp] = kde(radial_dist_last_full, 'Support', 'positive');
h(5) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [0 0.5 0], "FaceColor", [0 0.5 0]);
[fp,xfp] = kde(radial_dist_last_in', 'Support', 'positive');
h(4) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, 'LineStyle', ':', "EdgeColor", [0 0.8 0], "FaceColor", [0 0.8 0]);
[fp,xfp] = kde(radial_dist_first_full', 'Support', 'positive');
h(3) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, "EdgeColor", [0.5 0 0.5], "FaceColor", [0.5 0 0.5]);
[fp,xfp] = kde(mapped_radial_position_out, 'Support', 'positive');
h(2) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, 'LineStyle', '--', "EdgeColor", [1 0 1], "FaceColor", [1 0 1]);
[fp,xfp] = kde(radial_dist_last_out, 'Support', 'positive', 'Bandwidth', 7);
h(1) = area(xfp, fp, 'FaceAlpha', 0.2, "LineWidth", 2, 'LineStyle', '--', "EdgeColor", [0 1 0], "FaceColor", [0 1 0]);
ax = gca;
set(ax, 'FontSize', 34.5);
le = legend([h(3) h(5) h(6) h(4) h(2) h(1)], ...
    {sprintf('%.3gh (n=%d)', 36+10*(first_frame - 1)/60, length(radial_dist_first_full)), ...
    sprintf('%.3gh (n=%d)', 36+10*(last_frame - 1)/60, length(radial_dist_last_full)), ...
    sprintf('in %.3gh (n=%d)', 36+10*(first_frame - 1)/60, length(mapped_radial_position_in)), ...
    sprintf('in %.3gh (n=%d)', 36+10*(last_frame - 1)/60, length(radial_dist_last_in)), ...
    sprintf('out %.3gh (n=%d)', 36+10*(first_frame - 1)/60, length(mapped_radial_position_out)), ...
    sprintf('out %.3gh (n=%d)', 36+10*(last_frame - 1)/60, length(radial_dist_last_out))}, ...
    'Location','northwest', 'FontSize', 33, 'Box','off');
set(ax, "FontWeight", 'bold');
fontname(ax, "Arial");
axis square
set(ax, 'Units', 'centimeters', 'OuterPosition', [0 0 28 28]);
set(ax, "LineWidth", 4)
xlabel(ax, "r (µm)")
yticks([0, 0.04])
ylab = ylabel(ax, "density");
ylab.Position(1) = ylab.Position(1) + 40;
xticks([0, 100, 200, 300, 400]);
box off
exportgraphics(ax,'mapping_radial_distribution_in_out.png','Resolution',300)
close(f)
%% Angle of deviation for both inward and outward moving cells
angle_motion_in = [];
angle_motion_out = [];
for i = tracked_cell_ids
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    parent = cell_position_tracking(i).parent;
    daughter_last_pos = positions(end, :);
    pixel_pos = round(daughter_last_pos);
    parent_angle_motion = [];
    while parent ~= 0
        positions = cell_position_tracking(parent).positions;
        frames = cell_position_tracking(parent).frames;
        frame_first = frames(1);
        parent_angle_motion_p = [];
        if height(frames) > 1
            for k= 1:(height(frames) - 1)
                v1 = positions(k+1, :) - positions(k, :);
                v2 = colony_center - positions(1, :);
                angle = acos(dot(v1, v2)/(norm(v1)*norm(v2)));
                parent_angle_motion_p(end + 1, 1) = sign(det([v1; v2])) * angle;
            end
            parent_angle_motion = cat(1, parent_angle_motion, parent_angle_motion_p);
        end
        parent = cell_position_tracking(parent).parent;
    end   
    frame_earliest = frames(1, :);
    if frame_earliest > first_frame
        continue
    end
    positions = cell_position_tracking(i).positions;
    frames = cell_position_tracking(i).frames;
    if height(frames) > 1
        for k= 1:(height(frames) - 1)
            v1 = positions(k+1, :) - positions(k, :);
            v2 = colony_center - positions(1, :);
            angle = acos(dot(v1, v2)/(norm(v1)*norm(v2)));
            if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
                angle_motion_in(end + 1, 1) = sign(det([v1; v2])) * angle;
            end
            if mask_out_only(pixel_pos(2), pixel_pos(1)) ~= 0
                angle_motion_out(end + 1, 1) = sign(det([v1; v2])) * angle;
            end
        end
    end
    if mask_out_only(pixel_pos(2), pixel_pos(1)) == 0
        angle_motion_in = cat(1, angle_motion_in, parent_angle_motion);
    else
        angle_motion_out = cat(1, angle_motion_out, parent_angle_motion);
    end
end
%% Polarhistogram of angle of deviation for inward moving cells only
n = length(angle_motion_in);
bw = 2 * iqr(angle_motion_in) / n^(1/3);
n_bins = round((max(angle_motion_in) - min(angle_motion_in)) / bw);
f = figure;
f.Position = [10, 10, 2000, 2000];
pax = polaraxes;
polarhistogram(pax, angle_motion_in, n_bins + 2, "DisplayStyle", "bar", "FaceColor",[0, 1, 0], "LineWidth", 2);
pax.ThetaZeroLocation = 'top';
pax.ThetaLim = [-180, 180];
set(pax, "LineWidth", 4, "FontSize", 34.5)
fontname(pax, "Arial");
set(pax, 'Units', 'centimeters', 'OuterPosition', [0 0 26 26]);
rule = pax.RAxis;
rule.TickValues = 0:100:700;
set(pax, 'FontSize', 34.5, 'FontWeight', 'bold');
rule.Limits = [0, 650];
set(rule, "FontSize", 25);
rlabel = rule.Label;
rlabel.String = 'count';
rlabel.Rotation = 0;
rlabel.Position(2) = rlabel.Position(2) + 360;
rlabel.Position(1) = rlabel.Position(1) + 40;
thetaax = pax.ThetaAxis;
thetaax = thetaax.Label;
thetaax.String = 'angle of deviation';
pax.RTickLabel = []; 
for i = 0:100:600
    text(0, i, sprintf("%d", i), 'FontSize', 25, 'FontWeight', 'bold', 'FontName', 'Arial');
end
exportgraphics(pax, "polarhistogram_inward_only.png", 'Resolution', 300);
close(f)
%% Polarhistogram of angle of deviation for both inward and outward moving cells
n = length(angle_motion_in);
bw = 2 * iqr(angle_motion_in) / n^(1/3);
n_bins = round((max(angle_motion_in) - min(angle_motion_in)) / bw);
f = figure;
f.Position = [10, 10, 2000, 2000];
pax = polaraxes;
polarhistogram(pax, angle_motion_in, n_bins, "DisplayStyle", "stairs", "LineWidth", 2, 'Normalization', 'pdf');
hold on
n = length(angle_motion_out);
bw = 2 * iqr(angle_motion_out) / n^(1/3);
n_bins = round((max(angle_motion_out) - min(angle_motion_out)) / bw);
polarhistogram(pax, angle_motion_out, n_bins - 3, "DisplayStyle", "stairs", "LineWidth", 2, 'Normalization', 'pdf');
le = legend('in', 'out', 'Location','northeast', 'FontSize', 34.5, 'Box','off');
le.Position(1) = 0.45;
pax.ThetaZeroLocation = 'top';
pax.ThetaLim = [-180, 180];
set(pax, "LineWidth", 4, "FontSize", 34.5)
fontname(pax, "Arial");
set(pax, 'Units', 'centimeters', 'OuterPosition', [0 0 26 26]);
rule = pax.RAxis;
rule.TickValues = 0:0.05:0.25;
set(pax, 'FontSize', 34.5, 'FontWeight', 'bold');
set(rule, "FontSize", 23);
rlabel = rule.Label;
rlabel.String = 'density';
rlabel.Rotation = 0;
rlabel.Position(2) = rlabel.Position(2) + 0.14;
rlabel.Position(1) = rlabel.Position(1) + 43;
thetaax = pax.ThetaAxis;
thetaax = thetaax.Label;
thetaax.String = 'angle of deviation';
pax.RTickLabel = []; 
for i = 0:0.25:0.25
    text(0, i, sprintf("%g", i), 'FontSize', 25, 'FontWeight', 'bold', 'FontName', 'Arial');
end
exportgraphics(pax, "polarhistogram_in_out.png", 'Resolution', 300);
close(f)