rng(7)
% Settings
mice = {'A','D', 'E', 'F', 'G'}; %Mouse B removed
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding';
Fs = 15;
cutoff = 1.5;
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% Analysis parameters
food_types = {'Cheese', 'Carrot'};
pre_bout_frames = Fs * 4;
n_shuffles = 1000;

% Storage
modCellIdx = struct();  % Pos/neg/non indices
activityChanges = struct();  % Nx4 matrix per group

for f = 1:length(food_types)
    food = food_types{f};
    modCellIdx.(food).pos = [];
    modCellIdx.(food).neg = [];
    modCellIdx.(food).non = [];
    activityChanges.(food).pos = [];
    activityChanges.(food).neg = [];
end

for m = 1:length(mice)
    letter = mice{m};
    file_path = fullfile(base_path, ['Mouse_' letter], 'Week 2', 'Session 2', [letter '_hungry.mat']);

    if ~isfile(file_path)
        fprintf('Missing file for Mouse_%s. Skipping...\n', letter);
        continue;
    end

    fprintf('Processing Mouse_%s...\n', letter);
    data = load(file_path);

    % Filter data
    data_hunger = data.sumC_hunger;
    data_hunger_filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        data_hunger_filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end

    for f = 1:length(food_types)
        food = food_types{f};
        try
            consumption = getfield(data.timestamps, ['Mouse_' letter], 'Week2', 'Session2', food, 'Consumption');
        catch
            continue;
        end

        if isempty(consumption) || any(consumption == 0), continue; end

        start_frame = consumption(1);
        end_frame = consumption(2);
        if start_frame <= pre_bout_frames || end_frame > size(data_hunger_filtered, 2), continue; end

        bout_len = end_frame - start_frame + 1;
        shortened_len = floor(bout_len * 1);

        real_diffs = mean(data_hunger_filtered(:, start_frame:start_frame + shortened_len - 1), 2) - ...
                     mean(data_hunger_filtered(:, start_frame - pre_bout_frames:start_frame - 1), 2);

       % Shuffled null using circshift by 300*s
null_diffs = zeros(size(data_hunger_filtered, 1), n_shuffles);
for s = 1:n_shuffles
    shift_amt = 300 * s;
    shifted = circshift(data_hunger_filtered, [0, shift_amt]);

    null_eat = shifted(:, start_frame:start_frame + shortened_len - 1);
    null_pre = shifted(:, start_frame - pre_bout_frames:start_frame - 1);

    null_diffs(:, s) = mean(null_eat, 2) - mean(null_pre, 2);
end

        % Classify cells
        for c = 1:size(data_hunger_filtered, 1)
            hi = prctile(null_diffs(c,:), 85);
            lo = prctile(null_diffs(c,:), 15);
            if real_diffs(c) > hi
                modCellIdx.(food).pos(end+1) = c + (m-1)*size(data_hunger_filtered,1);
                activityChanges.(food).pos(end+1,1) = real_diffs(c);
            elseif real_diffs(c) < lo
                modCellIdx.(food).neg(end+1) = c + (m-1)*size(data_hunger_filtered,1);
                activityChanges.(food).neg(end+1,1) = real_diffs(c);
            else
                modCellIdx.(food).non(end+1) = c + (m-1)*size(data_hunger_filtered,1);
            end
        end
    end
end



% Setup
mice = {'A','D', 'E', 'F', 'G'};
foods_high = {'Cheese'};
foods_low = {'Carrot'};
% Set total cell count (for index offset)
n_cells = 100;

% Initialize pooled Δactivity values across all mice (cell-level)
high_vals_neg = [];
low_vals_neg = [];
high_vals_pos = [];
low_vals_pos = [];

% Loop over each mouse to collect non-redundant Δactivity for each group
for m = 1:length(mice)
    letter = mice{m};
    base = (m-1)*n_cells;
    
    % --- Negatively modulated cells ---
    high_cells = [];
    for f = 1:length(foods_high)
        food = foods_high{f};
        if isfield(modCellIdx.(food), 'neg')
            idx = modCellIdx.(food).neg;
            vals = activityChanges.(food).neg;
            % Get indices in this mouse's range
            mouseMask = idx > base & idx <= base + n_cells;
            these_cells = idx(mouseMask);
            these_vals = vals(mouseMask);
            % Avoid duplicates (across high-value foods)
            [unique_cells, ia] = unique(these_cells);
            high_cells = [high_cells; [unique_cells', these_vals(ia)]];
        end
    end
    high_vals_neg = [high_vals_neg; high_cells(:,2)];

    low_cells = [];
    for f = 1:length(foods_low)
        food = foods_low{f};
        if isfield(modCellIdx.(food), 'neg')
            idx = modCellIdx.(food).neg;
            vals = activityChanges.(food).neg;
            mouseMask = idx > base & idx <= base + n_cells;
            these_cells = idx(mouseMask);
            these_vals = vals(mouseMask);
            [unique_cells, ia] = unique(these_cells);
            low_cells = [low_cells; [unique_cells', these_vals(ia)]];
        end
    end
    low_vals_neg = [low_vals_neg; low_cells(:,2)];

    % --- Positively modulated cells ---
    high_cells_pos = [];
    for f = 1:length(foods_high)
        food = foods_high{f};
        if isfield(modCellIdx.(food), 'pos')
            idx = modCellIdx.(food).pos;
            vals = activityChanges.(food).pos;
            mouseMask = idx > base & idx <= base + n_cells;
            these_cells = idx(mouseMask);
            these_vals = vals(mouseMask);
            [unique_cells, ia] = unique(these_cells);
            high_cells_pos = [high_cells_pos; [unique_cells', these_vals(ia)]];
        end
    end
    high_vals_pos = [high_vals_pos; high_cells_pos(:,2)];

    low_cells_pos = [];
    for f = 1:length(foods_low)
        food = foods_low{f};
        if isfield(modCellIdx.(food), 'pos')
            idx = modCellIdx.(food).pos;
            vals = activityChanges.(food).pos;
            mouseMask = idx > base & idx <= base + n_cells;
            these_cells = idx(mouseMask);
            these_vals = vals(mouseMask);
            [unique_cells, ia] = unique(these_cells);
            low_cells_pos = [low_cells_pos; [unique_cells', these_vals(ia)]];
        end
    end
    low_vals_pos = [low_vals_pos; low_cells_pos(:,2)];
end

% --- Stats: NEGATIVELY modulated ---
[p_neg, h_neg] = ranksum(high_vals_neg, low_vals_neg);
fprintf('Ranksum test (NEG): p = %.4f\n', p_neg);

% --- Stats: POSITIVELY modulated ---
[p_pos, h_pos] = ranksum(high_vals_pos, low_vals_pos);
fprintf('Ranksum test (POS): p = %.4f\n', p_pos);

% --- Plot NEG ---
group_means = [nanmean(high_vals_neg), nanmean(low_vals_neg)];
group_sems = [nanstd(high_vals_neg)/sqrt(length(high_vals_neg)), ...
              nanstd(low_vals_neg)/sqrt(length(low_vals_neg))];

figure;
hold on;
c = bar(0.15, group_means(1), 'FaceColor', 'flat');
b = bar(0.35, group_means(2), 'FaceColor', 'flat');
c.CData = [1.0, 0.898, 0.6]; b.CData = [0.714, 0.894, 0.757];
c.EdgeColor = 'none'; b.EdgeColor = 'none';
c.BarWidth = 0.15; b.BarWidth = 0.15;
errorbar([0.15 0.35], group_means, group_sems, 'k.', 'LineWidth', 1);
%plot(0.35*ones(size(low_vals_neg)), low_vals_neg, 'ok', 'MarkerFaceColor', 'k');
%plot(0.15*ones(size(high_vals_neg)), high_vals_neg, 'ok', 'MarkerFaceColor', 'k');
set(gca,'linewidth',2, 'TickDir','out', 'xtick', [], 'fontsize', 24);
xlim([0 0.5]); ylim([-1.7 0]); yticks([-1.7 0]);
ylabel({'mean df/F (z-scored)'; '(post-pre)'}, 'FontSize', 24);
legend({'high value foods', 'low value foods'}, 'Location', 'eastoutside'); legend boxoff; box off;

% Significance stars
aster = "";
if p_neg < 0.05 && p_neg >= 0.01, aster = "*";
elseif p_neg < 0.01 && p_neg >= 0.001, aster = "**";
elseif p_neg < 0.001 && p_neg > 0.0001, aster = "***";
elseif p_neg < 0.0001, aster = "****";
end
text(0.24, 0.05, aster, 'HorizontalAlignment', 'center', 'FontSize', 30);
text(0.24, 0.1, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

% --- Plot POS ---
group_means_pos = [nanmean(high_vals_pos), nanmean(low_vals_pos)];
group_sems_pos = [nanstd(high_vals_pos)/sqrt(length(high_vals_pos)), ...
                  nanstd(low_vals_pos)/sqrt(length(low_vals_pos))];

figure;
hold on;
c = bar(0.15, group_means_pos(1), 'FaceColor', 'flat');
b = bar(0.35, group_means_pos(2), 'FaceColor', 'flat');
c.CData = [1.0, 0.898, 0.6]; b.CData = [0.714, 0.894, 0.757];
c.EdgeColor = 'none'; b.EdgeColor = 'none';
c.BarWidth = 0.15; b.BarWidth = 0.15;
errorbar([0.15 0.35], group_means_pos, group_sems_pos, 'k.', 'LineWidth', 1);
%plot(0.35*ones(size(low_vals_pos)), low_vals_pos, 'ok', 'MarkerFaceColor', 'k');
%plot(0.15*ones(size(high_vals_pos)), high_vals_pos, 'ok', 'MarkerFaceColor', 'k');
set(gca,'linewidth',2, 'TickDir','out', 'xtick', [], 'fontsize', 24);
xlim([0 0.5]); ylim([0 1.7]); yticks([0 1.7]);
ylabel({'mean df/F (z-scored)'; '(post-pre)'}, 'FontSize', 24);
legend({'high value foods', 'low value foods'}, 'Location', 'eastoutside'); legend boxoff; box off;

% Significance stars
aster = "";
if p_pos < 0.05 && p_pos >= 0.01, aster = "*";
elseif p_pos < 0.01 && p_pos >= 0.001, aster = "**";
elseif p_pos < 0.001 && p_pos > 0.0001, aster = "***";
elseif p_pos < 0.0001, aster = "****";
end
text(0.24, 1.65, aster, 'HorizontalAlignment', 'center', 'FontSize', 30);
%text(0.24, 1.70, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);







%%%%%%%%%%%%%%%%%%%%%%%%% 
% Example heatmaps from Mouse B, Week 2, Session 2 (this first block of
% code is heatmap of all cells
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_B/Week 2/Session 2/B_hungry.mat");

% --- Parameters ---
startFrame = timestamps.Mouse_B.Week2.Session2.Cheese.Consumption(1); %This makes it so that the heatmap is sorted from the first frame of Cheese Bout #1
endFrame = timestamps.Mouse_B.Week2.Session2.Cheese.Consumption(2);  % Or specify an upper bound like 12000
windowSize = 30;  % e.g., 30 frames = 2 seconds at 15 Hz

% --- Extract relevant window ---
segment = sumC_hunger(:, startFrame:endFrame);

% --- Preallocate vector for max average per cell ---
numCells = size(segment, 1);
maxAvgActivity = zeros(numCells, 1);

% --- For each cell, compute max average across sliding windows ---
for i = 1:numCells
    trace = segment(i, :);
    slidingAvg = movmean(trace, windowSize, 'Endpoints', 'discard');
    maxAvgActivity(i) = max(slidingAvg);
end

% --- Sort by peak average activity ---
[~, sort_order] = sort(maxAvgActivity, 'descend');  % or 'ascend' if preferred

% --- Apply sorting to full matrix ---
sorted_sumC = sumC_hunger(sort_order, :);

% --- Plot heatmap ---
figure;
imagesc(sorted_sumC);
colormap parula;
colorbar;
clim([-1.5 3]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Crop parameters --- THIS IS MOUSE B, WEEK 2, SESSION 2, negative
% cells
crop_start = 10000;
crop_end = 25000;
fps = 15;

% --- Extract cropped data ---
cropped_data = sorted_sumC(40:end, crop_start:crop_end);
num_cells = size(cropped_data, 1);
num_frames = size(cropped_data, 2);

% --- Plot heatmap ---
figure;
imagesc(cropped_data);
set(gca,'YDir', 'normal')
colormap parula;
c= colorbar;

hold on;
clim([-1.5 3]);
% Get current color limits
c_min = clim;  % [min, max]
% Set tick positions at the min and max
c.Ticks = [c_min(1), c_min(2)];

% Set labels to 'min' and 'max'
c.TickLabels = {'min', 'max'};

% Optional styling
c.FontSize = 20;
% --- X-axis: convert frames to minutes ---
xticks_in_frames = linspace(1, num_frames, 6);
xticklabels_in_min = round(((xticks_in_frames + crop_start - 1) / fps) / 60, 1);

xlabel('time (min)', 'FontSize', 30);

% --- Y-axis formatting ---
yticks([1 num_cells]);
yticklabels({'1', num2str(num_cells)});
ylabel('cell #', 'FontSize', 30);

% --- Cheese periods (yellow-beige bar) ---
cheese_times = timestamps.Mouse_B.Week2.Session2.Cheese.Consumption;
for i = 1:2:length(cheese_times)
    x1 = cheese_times(i) - crop_start + 1;
    x2 = cheese_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -1, x2 - x1, 1], ...
                  'FaceColor', [1.0, 0.898, 0.6], ...
                  'EdgeColor', 'none');
    end
end

% --- Broccoli periods (green bar) ---
broccoli_times = timestamps.Mouse_B.Week2.Session2.Broccoli.Consumption;
for i = 1:2:length(broccoli_times)
    x1 = broccoli_times(i) - crop_start + 1;
    x2 = broccoli_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -1, x2 - x1, 1], ...
                  'FaceColor', [0.714, 0.894, 0.757], ...
                  'EdgeColor', 'none');
    end
end
% === Remove x-axis labels ===
xticks([]);  % removes x-ticks
xlabel('');  % removes x-axis label

% === Add white scale bar (900 frames, 60 sec) ===
bar_length = 900;  % frame
bar_height = 0.5;    % vertical thickness of the bar in cell units

% Position: bottom right
x_start = size(cropped_data, 2) - bar_length - 200;  % 20-frame margin
y_start = num_cells - 2.5;  % just below the heatmap

% Draw white rectangle as scale bar
rectangle('Position', [x_start, y_start-0.5, bar_length, bar_height], ...
          'FaceColor', 'w', ...
          'EdgeColor', 'none');

hold on;
% Add the text label centered below the bar
text(x_start+150, y_start, '60 s', ...
     'Color', 'w', ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 20, ...
     'FontWeight', 'bold');
% --- Adjust axes ---
ylim([-1.5 num_cells]);
set(gca, 'FontSize', 30);
box off;
exportgraphics(gca,'cheeseVSbroc_negative.eps','ContentType','vector')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Crop parameters --- THIS IS MOUSE B, WEEK 2, SESSION 2, positive
% cells
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_A/Week 2/Session 2/A_hungry.mat");
%%%%%%%%%%%%%%%%%%%%%%%% 
%Example heatmaps from Mouse B, Week 2, Session 2 (this first block of
%code is heatmap of all cells
%--- Parameters ---
startFrame = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption(1); %This makes it so that the heatmap is sorted from the first frame of Cheese Bout #1
endFrame = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption(2);  % Or specify an upper bound like 12000
windowSize = 30;  % e.g., 30 frames = 2 seconds at 15 Hz

%--- Extract relevant window ---
segment = sumC_hunger(:, startFrame:endFrame);

%--- Preallocate vector for max average per cell ---
numCells = size(segment, 1);
maxAvgActivity = zeros(numCells, 1);

%--- For each cell, compute max average across sliding windows ---
for i = 1:numCells
    trace = segment(i, :);
    slidingAvg = movmean(trace, windowSize, 'Endpoints', 'discard');
    maxAvgActivity(i) = max(slidingAvg);
end

%--- Sort by peak average activity ---
[~, sort_order] = sort(maxAvgActivity, 'descend');  % or 'ascend' if preferred

%--- Apply sorting to full matrix ---
sorted_sumC = sumC_hunger(sort_order, :);
sorted_sumC(9:10,:) = [];
sorted_sumC(5,:) = [];
sorted_sumC(1,:) = [];
sorted_sumC(10,:) = [];

crop_start = 10000;
crop_end = 25000;
fps = 15;
cropped_data = sorted_sumC(1:14, crop_start:crop_end);
num_cells = size(cropped_data, 1);
num_frames = size(cropped_data, 2);

%--- Plot heatmap ---
figure;
imagesc(cropped_data);
colormap parula;
c = colorbar;
clim([-1.5 3]);
hold on;
% Get current color limits
c_min = clim;  % [min, max]
% Set tick positions at the min and max
c.Ticks = [c_min(1), c_min(2)];

% Set labels to 'min' and 'max'
c.TickLabels = {'min', 'max'};

% Optional styling
c.FontSize = 20;
% --- X-axis: convert frames to minutes ---
xticks_in_frames = linspace(1, num_frames, 6);
xticklabels_in_min = round(((xticks_in_frames + crop_start - 1) / fps) / 60, 1);

xlabel('time (min)', 'FontSize', 30);

% --- Y-axis formatting ---
yticks([1 num_cells]);
yticklabels({'1', num2str(num_cells)});
ylabel('cell #', 'FontSize', 30);

% --- Cheese periods (yellow-beige bar) ---
cheese_times = timestamps.Mouse_A.Week2.Session2.Cheese.Consumption;
for i = 1:2:length(cheese_times)
    x1 = cheese_times(i) - crop_start + 1;
    x2 = cheese_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -0.25, x2 - x1, 0.5], ...
                  'FaceColor', [1.0, 0.898, 0.6], ...
                  'EdgeColor', 'none');
    end
end

% --- Broccoli periods (green bar) ---
broccoli_times = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption;
for i = 1:2:length(broccoli_times)
    x1 = broccoli_times(i) - crop_start + 1;
    x2 = broccoli_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -0.25, x2 - x1, 0.5], ...
                  'FaceColor', [0.714, 0.894, 0.757], ...
                  'EdgeColor', 'none');
    end
end
% === Remove x-axis labels ===
xticks([]);  % removes x-ticks
xlabel('');  % removes x-axis label

% === Add white scale bar (900 frames, 60 sec) ===
bar_length = 900;  % frame
bar_height = 0.25;    % vertical thickness of the bar in cell units

% Position: bottom right
x_start = size(cropped_data, 2) - bar_length - 200;  % 20-frame margin
y_start = num_cells - 1.5;  % just below the heatmap

% Draw white rectangle as scale bar
rectangle('Position', [x_start, y_start+0.2, bar_length, bar_height], ...
          'FaceColor', 'w', ...
          'EdgeColor', 'none');

hold on;
% Add the text label centered below the bar
text(x_start+150, y_start+0.4, '60 s', ...
     'Color', 'w', ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 20, ...
     'FontWeight', 'bold');
% --- Adjust axes ---
ylim([-0.5 num_cells]);
set(gca, 'FontSize', 30);
box off;
exportgraphics(gca,'cheesevsbroc_positive.eps','ContentType','vector')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_A/Week 2/Session 2/A_hungry.mat");
%%%%%%%%%%%%%%%%%%%%%%%% 
%Example heatmaps from positively modulated cells from Mouse A, Week 2, Session 2
%--- Parameters ---
startFrame = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption(1); %This makes it so that the heatmap is sorted from the first frame of Cheese Bout #1
endFrame = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption(2);  % Or specify an upper bound like 12000
windowSize = 30;  % e.g., 30 frames = 2 seconds at 15 Hz

%--- Extract relevant window ---
segment = sumC_hunger(:, startFrame:endFrame);

%--- Preallocate vector for max average per cell ---
numCells = size(segment, 1);
maxAvgActivity = zeros(numCells, 1);

%--- For each cell, compute max average across sliding windows ---
for i = 1:numCells
    trace = segment(i, :);
    slidingAvg = movmean(trace, windowSize, 'Endpoints', 'discard');
    maxAvgActivity(i) = max(slidingAvg);
end

%--- Sort by peak average activity ---
[~, sort_order] = sort(maxAvgActivity, 'descend');  % or 'ascend' if preferred

%--- Apply sorting to full matrix ---
sorted_sumC = sumC_hunger(sort_order, :);
sorted_sumC(9:10,:) = [];
sorted_sumC(5,:) = [];
sorted_sumC(1,:) = [];
sorted_sumC(10,:) = [];

crop_start = 10000;
crop_end = 25000;
fps = 15;
cropped_data = sorted_sumC(1:14, crop_start:crop_end);
num_cells = size(cropped_data, 1);
num_frames = size(cropped_data, 2);

%--- Plot heatmap ---
figure;
imagesc(cropped_data);
colormap parula;
c = colorbar;
clim([-1.5 3]);
hold on;
% Get current color limits
c_min = clim;  % [min, max]
% Set tick positions at the min and max
c.Ticks = [c_min(1), c_min(2)];

% Set labels to 'min' and 'max'
c.TickLabels = {'min', 'max'};

% Optional styling
c.FontSize = 20;
% --- X-axis: convert frames to minutes ---
xticks_in_frames = linspace(1, num_frames, 6);
xticklabels_in_min = round(((xticks_in_frames + crop_start - 1) / fps) / 60, 1);

xlabel('time (min)', 'FontSize', 30);

% --- Y-axis formatting ---
yticks([1 num_cells]);
yticklabels({'1', num2str(num_cells)});
ylabel('cell #', 'FontSize', 30);

% --- Cheese periods (yellow-beige bar) ---
cheese_times = timestamps.Mouse_A.Week2.Session2.Cheese.Consumption;
for i = 1:2:length(cheese_times)
    x1 = cheese_times(i) - crop_start + 1;
    x2 = cheese_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -0.25, x2 - x1, 0.5], ...
                  'FaceColor', [1.0, 0.898, 0.6], ...
                  'EdgeColor', 'none');
    end
end

% --- Broccoli periods (green bar) ---
broccoli_times = timestamps.Mouse_A.Week2.Session2.Broccoli.Consumption;
for i = 1:2:length(broccoli_times)
    x1 = broccoli_times(i) - crop_start + 1;
    x2 = broccoli_times(i+1) - crop_start + 1;
    if x1 > 0 && x2 <= num_frames
        rectangle('Position', [x1, -0.25, x2 - x1, 0.5], ...
                  'FaceColor', [0.714, 0.894, 0.757], ...
                  'EdgeColor', 'none');
    end
end
% === Remove x-axis labels ===
xticks([]);  % removes x-ticks
xlabel('');  % removes x-axis label

% === Add white scale bar (900 frames, 60 sec) ===
bar_length = 900;  % frame
bar_height = 0.25;    % vertical thickness of the bar in cell units

% Position: bottom right
x_start = size(cropped_data, 2) - bar_length - 200;  % 20-frame margin
y_start = num_cells - 1.5;  % just below the heatmap

% Draw white rectangle as scale bar
rectangle('Position', [x_start, y_start+0.2, bar_length, bar_height], ...
          'FaceColor', 'w', ...
          'EdgeColor', 'none');

hold on;
% Add the text label centered below the bar
text(x_start+150, y_start+0.4, '60 s', ...
     'Color', 'w', ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', ...
     'FontSize', 20, ...
     'FontWeight', 'bold');
% --- Adjust axes ---
ylim([-0.5 num_cells]);
set(gca, 'FontSize', 30);
box off;
exportgraphics(gca,'cheesevsbroc_positive.eps','ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can manipulate this code to look at activity in different ranges
% Currently showing cheese/broccoli eating in Mouse A hungry
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_A/Week 2/Session 2/A_hungry.mat");

startFrame = timestamps.Mouse_A.Week2.Session2.Cheese.Consumption(1); %This makes it so that the heatmap is sorted from the first frame of Cheese Bout #1
endFrame = timestamps.Mouse_A.Week2.Session2.Cheese.Consumption(2);  % Or specify an upper bound like 12000
windowSize = 30;  % e.g., 30 frames = 2 seconds at 15 Hz

% --- Extract relevant window ---
segment = sumC_hunger(:, startFrame:endFrame);

% --- Preallocate vector for max average per cell ---
numCells = size(segment, 1);
maxAvgActivity = zeros(numCells, 1);

% --- For each cell, compute max average across sliding windows ---
for i = 1:numCells
    trace = segment(i, :);
    slidingAvg = movmean(trace, windowSize, 'Endpoints', 'discard');
    maxAvgActivity(i) = max(slidingAvg);
end

% --- Sort by peak average activity ---
[~, sort_order] = sort(maxAvgActivity, 'descend');  % or 'ascend' if preferred

% --- Apply sorting to full matrix ---
sorted_sumC = sumC_hunger(sort_order, :);



% --- Plot heatmap ---
figure;
imagesc(sorted_sumC(:, 10000:25000));
colormap parula;
colorbar;
clim([-1.5 3]);

Fs = 15;
cutoff = 1.5;
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

data_hunger = sorted_sumC;
    data_hunger_filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        data_hunger_filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end


% --- Parameters ---
fps = 15;
window_sec = 4;
window_frames = window_sec * fps;

% --- Data (remove first 29 cells as specified) ---
data = data_hunger_filtered(60:end, :);  % [cells x time]

% --- Center time point ---
center_frame = timestamps.Mouse_D.Week2.Session2.Cheese.Consumption(1);

% --- Time window around the center ---
start_frame = center_frame - window_frames;
end_frame = center_frame + window_frames;

% --- Check bounds ---
if start_frame < 1 || end_frame > size(data, 2)
    error('Window exceeds data bounds');
end

% --- Extract data window ---
activity_window = data(:, start_frame:end_frame);  % [cells x 2*window_frames+1]

% --- Compute average across cells ---
mean_activity = mean(activity_window, 1);
sem_activity = std(activity_window, 0, 1) / sqrt(size(activity_window, 1));

% --- Time axis in seconds (centered at 0) ---
t = linspace(-window_sec, window_sec, size(activity_window, 2));
figure;
% Plot SEM as a shaded region
fill([t fliplr(t)], ...
     [mean_activity + sem_activity, fliplr(mean_activity - sem_activity)], ...
     [1.0, 0.898, 0.6], ...       % light yellow
     'FaceAlpha', 0.3, ...
     'EdgeColor', 'none');
hold on;

% --- Plot ---

plot(t, mean_activity, 'Color',[1.0, 0.898, 0.6], 'LineWidth', 2);
xlim([-4 4]);
%ylim([-1.1 0.5]);
%yticks([-1.1 0.5]);
%yticklabels({"-1.1", "0.5"});
xticks([-4 0 4]);
xticklabels({"-4","0","4"});
xline(0, "r--", "LineWidth",2);
yline(0, "Color", [0.5,0.5,0.5],"LineWidth",2);
set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 30);
xlabel('time from eating onset (s)', 'FontSize', 30);
ylabel('mean df/f (z-scored)', 'FontSize', 30);
grid off;
box off;