% === Load data ===
data_dir = '/Users/nitishpatel/Documents/MATLAB/GCaMP_CheeseVsBroc';
load(fullfile(data_dir, 'eating_time.mat'));   % contains eating_time.<Mouse>.Cheese/Broccoli

% === Parameters ===
Fs = 7.5;
window_sec = 5;  % last 5 seconds
window_frames = round(window_sec * Fs);

% === Get all mouse names in structure ===
mouse_names = {'FP1','FP3','FP5','FP6','FP7','FP8'};
nMice = numel(mouse_names);

% === Initialize output arrays ===
auc_cheese_all = nan(nMice, 1);
auc_broc_all   = nan(nMice, 1);

% === Loop through each mouse ===
for i = 1:nMice
    mouse = mouse_names{i};
    fprintf('Processing %s...\n', mouse);
  
    % Load mouse-specific photometry file
    fp_file = fullfile(data_dir, [mouse '_000.mat']);
    if ~isfile(fp_file)
        warning('File not found: %s, skipping...', fp_file);
        continue;
    end
    load(fp_file, 'sig', 'ref');  % load only relevant variables

    % Preprocess photometry
    if strcmp(mouse, 'FP5')
       sig(7749:7750) = sig(7747:7748); %FP5
       sig = sig(900:end); %FP5, FP7
       ref = ref(900:end); %FP5, FP7
    end
    if strcmp(mouse, 'FP6')
       sig = sig(900:end); %FP6
       ref = ref(900:end); %FP6
    end
    sig(end) = sig(end - 1);
    ref(end) = ref(end - 1);

    p = polyfit(ref, sig, 1);
    s = sig - (ref * p(1) + p(2));

    s = zscore(s);

    % Extract bout end frame indices
    get_bout_end_indices = @(mat) round((mat(:,3)*60 + mat(:,4)) * Fs);

    try
        cheese_ends = get_bout_end_indices(eating_time.(mouse).Cheese);
        broccoli_ends = get_bout_end_indices(eating_time.(mouse).Broccoli);
    catch
        warning('Missing or malformed bout info for %s, skipping...', mouse);
        continue;
    end

    %% === MODIFIED: Equalize bout count across foods ===
    if isempty(cheese_ends) || isempty(broccoli_ends)
        continue;
    end
    n_use = min(length(cheese_ends), length(broccoli_ends));

    rng(0);  % for reproducibility
    cheese_ends = cheese_ends(randperm(length(cheese_ends), n_use));
    broccoli_ends = broccoli_ends(randperm(length(broccoli_ends), n_use));

    % Compute average AUC for last 5 seconds of each bout (equalized)
    compute_avg_auc = @(ends) mean(arrayfun(@(idx) ...
        trapz(s(max(1, idx - window_frames + 1):min(idx, length(s)))), ends));

    auc_cheese_all(i) = compute_avg_auc(cheese_ends);
    auc_broc_all(i)   = compute_avg_auc(broccoli_ends);
end

% === Remove NaNs (mice without data) ===
valid_idx = ~isnan(auc_cheese_all) & ~isnan(auc_broc_all);
auc_cheese = auc_cheese_all(valid_idx);
auc_broc = auc_broc_all(valid_idx);

% === Plotting ===
means = [mean(auc_cheese), mean(auc_broc)];
sems = [std(auc_cheese)/sqrt(length(auc_cheese)), ...
        std(auc_broc)/sqrt(length(auc_broc))];

figure;
hold on;

b = bar(0.15, [mean(auc_cheese)], 'FaceColor', 'flat');
hold on;
c = bar(0.35, [mean(auc_broc)], 'FaceColor', 'flat');
b.CData = [1.00, 0.92, 0.60];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
b.FaceAlpha = 0.6;
hold on;
c.CData = [0.75, 1, 0.75]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
c.FaceAlpha = 0.6;

% Add error bars
errorbar([0.15 0.35], means, sems, 'k.', 'LineWidth', 1);

% Add individual mouse dots
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
plot(0.15*ones(size(auc_cheese)), auc_cheese, 'ok', 'MarkerFaceColor', 'k');
plot(0.35*ones(size(auc_broc)), auc_broc, 'ok', 'MarkerFaceColor', 'k');
xlim([0 0.5])
ylim([-70 0])  
yticks([-70 0]);
ylabel('mean AUC of df/F (z-scored)');
legend({'cheese', 'broccoli'}, 'Location', 'eastoutside');
legend boxoff;
box off;

[p, h] = signrank(auc_cheese, auc_broc);
aster = "";
if p < 0.05 && p >= 0.01
    aster = "*";
elseif p < 0.01 && p >= 0.001
    aster = "**";
elseif p < 0.001 && p > 0.0001
    aster = "***";
elseif p < 0.0001
    aster = "****";
end

text(0.24, 3,  aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 5, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
%exportgraphics(gca,'FP_cheesevsbroc.eps','ContentType','vector')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below is example trace from FP5
% === Data ===
fps = 7.5;

% === Duration in frames ===
ten_min_frames = round(10 * 60 * fps);  % 10 minutes worth of frames

% === Frame indices for each 10-min segment ===
idx1 = 1:25:ten_min_frames;
idx2 = ten_min_frames+1:25:2*ten_min_frames;
idx3 = 2*ten_min_frames+1:25:max(length(s), 3*ten_min_frames);

% === Time vector for x-axis ===
t = (0:length(s)-1) / fps / 60;  % in minutes

% === Plot ===
figure; hold on;
plot(t(idx1), s(idx1), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);       % dark gray
plot(t(idx2), s(idx2), 'Color', [0 0.5 0], 'LineWidth', 1.5);           % dark green
plot(t(idx3), s(idx3), 'Color', [0.6 0.5 0], 'LineWidth', 1.5);         % dark yellow
plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);

% === Add transparent yellow bout rectangles ===
cheese_bout_times = [
    20 + 16/60, 21 + 37/60;  % Note that 18 seconds had to be added to scored timestamps (eating_time) to align with neural activity
    22 + 24/60, 23 + 36/60;
    24 + 12/60, 25 + 27/60;
    26 +  3/60, 27 + 50/60;
];

yl = ylim;
for i = 1:size(cheese_bout_times, 1)
    x1 = cheese_bout_times(i, 1);
    x2 = cheese_bout_times(i, 2);
    rectangle('Position', [x1, yl(1), x2 - x1, yl(2) - yl(1)], ...
              'FaceColor', [1, 0.92, 0.6, 0.3], 'EdgeColor', 'none');
end


broc_bout_times = [
    11 + 4/60, 12 + 22/60;  % Note that 18 seconds had to be added to scored timestamps (eating_time) to align with neural activity
    12 + 56/60, 13 + 7/60;
    13 + 25/60, 13 + 35/60;
    13 + 40/60, 14 + 11/60;
    14 +  57/60, 15 + 12/60;
];

yl = ylim;
for i = 1:size(broc_bout_times, 1)
    x1 = broc_bout_times(i, 1);
    x2 = broc_bout_times(i, 2);
    rectangle('Position', [x1, yl(1), x2 - x1, yl(2) - yl(1)], ...
              'FaceColor', [0.75, 1, 0.75, 0.3], 'EdgeColor', 'none');
end

xlabel('time (min)');
xlim([0 30])
xticks([0 10 20 30])
xticklabels({'0', '10', '20', '30'})
yticks([-3 0 3]);
yticklabels([-3 0 3]);
ylabel('example df/F (z-scored)');
set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 24)
% === Create dummy plots for legend ===
h_baseline = plot(nan, nan, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);     % baseline line
h_broc     = patch(nan, nan, [0.75, 1, 0.75], 'FaceAlpha', 0.6, 'EdgeColor', 'none'); % pastel green box
h_cheese   = patch(nan, nan, [1, 0.92, 0.6], 'FaceAlpha', 0.6, 'EdgeColor', 'none');  % pastel yellow box

legend([h_baseline, h_broc, h_cheese], {'baseline', 'broccoli eating', 'cheese eating'}, ...
    'Location', 'eastoutside');
legend boxoff;
legend boxoff;
box off;
%exportgraphics(gca,'FP5_cheesevsbroc.eps','ContentType','vector')
