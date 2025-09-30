% Define paths
%To run this code: you need [Letter]_ghrelin_coreg.mat file with
%injection_frame_ghrelin and injection_frame_saline as the frame number at
%which injection occurred (can find this using the code above)


%IMPORTANT: For correlation graph: use "_ghrelin_coreg.mat",
%"_ghrelin.mat", and "_saline.mat". AND use ghrelin_pre: 1:9000,
%ghrelin_post: size(data,2) - 8100:size(data,2);
% use saline_pre: 2000:9000; use saline_post: size(data_sal, 2) - 8100:size(data_sal, 2)

%IMPORTANT: For bar graph: use "_ghrelin_coreg(new2).mat", "_ghrelin.mat",
%and "_saline(old).mat". AND use ghrelin_pre: 1:10000, ghrelin_post:
%10000:15000; use saline_pre: 1900:9800, and saline_post: 51000:59000


base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin';
mouse_list = {'Mouse_A', 'Mouse_F'};
letters = {'A', 'F'};
injection_frame_ghrelin = 10000;
injection_frame_saline = 9500;
% Initialize result storage
all_post_minus_pre = [];

for i = 1:length(mouse_list)
    mouse_folder = fullfile(base_path, mouse_list{i});
    letter = letters{i};
    coreg_file = fullfile(mouse_folder, [letter '_ghrelin_coreg.mat']);
    
    if isfile(coreg_file)
        fprintf('Processing %s...\n', letter);
        
        % Load co-registration file
        load(coreg_file, 'cell_registered_struct');
        A = cell_registered_struct.cell_to_index_map;
        A(any(A == 0, 2), :) = [];
        rows = A(:,1);
        rows_sal = A(:,2);
        
        % Load ghrelin and saline trace data
        ghrelin_file = fullfile(mouse_folder, [letter '_ghrelin.mat']);
        saline_file = fullfile(mouse_folder, [letter '_saline.mat']); 
        
        if ~isfile(ghrelin_file) || ~isfile(saline_file)
            warning('Missing ghrelin or saline .mat file for %s. Skipping.', letter);
            continue;
        end
        
        load(ghrelin_file, 'sumC'); % ghrelin data
        data = sumC(rows, :);
        %data(11:12,:) = [];
        load(saline_file, 'sumC_sal'); % saline data
        data_sal = sumC_sal(rows_sal, :);
        %data_sal(11:12,:) = [];
        n_cells = size(data, 1);
        post_minus_pre = nan(n_cells, 2); % column 1 = ghrelin, column 2 = saline

         % ----------------
        % GHRELIN ANALYSIS
        % ----------------
        window_ghrelin = 9000; % 9 min window
        pre_start_ghrelin = 1;
        pre_end_ghrelin = 9000; %10000;
       %post_length_ghrelin = %size(data, 2) - injection_frame_ghrelin;
        post_start_ghrelin = size(data,2) - 8100; %10000
        post_end_ghrelin = size(data,2); %15000

        if post_end_ghrelin <= size(data,2)
            for c = 1:n_cells
                pre_mean = mean(data(c, pre_start_ghrelin:pre_end_ghrelin));
                post_mean = mean(data(c, post_start_ghrelin:post_end_ghrelin));
                post_minus_pre(c, 1) = post_mean - pre_mean;
            end
        else
            warning('Ghrelin window exceeds data length for %s. Skipping ghrelin.', letter);
        end

        % ----------------
        % SALINE ANALYSIS
        % ----------------
        window_saline = 7650; %8.5 min window
        pre_start_saline = 2000;%1900
        pre_end_saline = 9000;%9800
        %post_length_saline = size(data_sal, 2) - injection_frame_saline;
        post_start_saline = size(data_sal, 2) - 8100; %51000
        post_end_saline = size(data_sal, 2); %59000

        if post_end_saline <= size(data_sal,2)
            for c = 1:n_cells
                pre_mean = mean(data_sal(c, pre_start_saline:pre_end_saline));
                post_mean = mean(data_sal(c, post_start_saline:post_end_saline));
                post_minus_pre(c, 2) = post_mean - pre_mean;
            end
        else
            warning('Saline window exceeds data length for %s. Skipping saline.', letter);
        end

        % Append to full matrix
        all_post_minus_pre = [all_post_minus_pre; post_minus_pre];
    else
        fprintf('No co-registration file for %s. Skipping.\n', letter);
    end
end

% all_post_minus_pre now contains the injection modulation changes for all mice
% Column 1 = ghrelin modulation; Column 2 = saline modulation



% Remove any rows with NaNs before plotting
valid_idx = all(~isnan(all_post_minus_pre), 2);
data_clean = all_post_minus_pre(valid_idx, :);

% Extract ghrelin and saline modulation
ghrelin_mod = data_clean(:, 1);
saline_mod = data_clean(:, 2);
%ghrelin_mod(15) = [];
%saline_mod(15) = [];

% Compute correlation
[r, p] = corr(ghrelin_mod, saline_mod);

% Create scatter plot
figure;
scatter(ghrelin_mod, saline_mod, 90, 'filled', "MarkerFaceColor", [0.63, 0.47, 0.35],"LineWidth",1.2, "MarkerEdgeAlpha", 0.8);
xlabel({'ghrelin df/F (z-scored)'; '(post-pre)'});
ylabel({'saline df/F (z-scored)'; '(post-pre)'});
text(-1.5, 0.8, sprintf('r = %.2f\np = %.1e', r, p), ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'center');

pfit = polyfit(ghrelin_mod, saline_mod, 1);  % Returns [slope, intercept]
% Define x values using the axis limits
x_vals = xlim;  % or [min(ghrelin_mod), max(ghrelin_mod)] if you prefer
y_vals = polyval(pfit, x_vals);  % Evaluate the line at these x values
hold on;
% Plot the custom regression line
plot(x_vals, y_vals, '-', 'Color', [0.63, 0.47, 0.35], 'LineWidth', 1.5);

set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 24)
xlim([-2.5 1.5])
xticks([-2.5 1.5])
yticks([-2.5 1.5])
ylim([-2.5 1.5])


grid off;
box off;
%exportgraphics(gca,'ghrelin_saline_corr.eps','ContentType','vector')


% Ensure data is clean (no NaNs)
valid_idx = all(~isnan(all_post_minus_pre), 2);
data_clean = all_post_minus_pre(valid_idx, :);
ghrelin_mod = data_clean(:, 1);
saline_mod = data_clean(:, 2);

figure;
b = bar(0.35, [mean(ghrelin_mod)], 'FaceColor', 'flat');
hold on;
c = bar(0.15, [mean(saline_mod)], 'FaceColor', 'flat');
b.CData = [65/255 111/255 178/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [249/255 218/255 120/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.15 0.35], [mean(saline_mod), mean(ghrelin_mod)], [std(saline_mod)/sqrt(length(saline_mod)), std(ghrelin_mod)/sqrt(length(ghrelin_mod))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
%plot(0.15, ch_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
%plot(0.35, yfp_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
xlim([0 0.5])
yticks([-0.8 0.4])
ylim([-0.8 0.4])
legend({'ghrelin', 'saline'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'mean df/F (z-scored)'; '(post-pre)'}, 'FontSize', 24);
%exportgraphics(gca,'MS_Ghrelin_Bar.eps','ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below is for plotting Ghrelin cell #9 and Saline cell #12 from Mouse_A as example traces
figure;
subplot(2, 1, 1);
avg_activity = data(9,:);
avg_activity = avg_activity(:, 1:54108);
frame_rate = 15;  % frames per second

% Convert x-values (frame indices) to time in minutes
time_minutes = (1:80:length(avg_activity)) / frame_rate / 60;

% Plot
plot(time_minutes, avg_activity(1:80:end), '-', 'LineWidth', 2, 'Color',[65/255 111/255 178/255]);
hold on;

% Add gray horizontal 0 line
plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on;
xline(10.5789, 'r--', 'LineWidth', 2);

% Darken axes
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% Set limits and ticks
xlim([0, 60]);  % ~0 to 5 seconds in minutes
ylim([-2, 7]);
xticks([0,30, 60]);  % 0s, 2.5s, 5s in minutes
%xticklabels({'0', '2.5', '5'});  % Optional: keep in seconds for clarity

yticks([-2,0, 7]);
yticklabels({'-2', '0', '7'});

% Aesthetics
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);

xlabel('time (minutes)');  % or 'Time (s)' if you keep ticks in seconds
ylabel('mean df/F (z-scored)');
legend off;
box off;
hold off;

subplot(2, 1, 2);
avg_saline = data_sal(12,:);
time_saline = (1:80:length(avg_saline)) / frame_rate / 60;
plot(time_saline, avg_saline(1:80:end), '-', 'LineWidth', 2, 'Color',[231/255 173/255 82/255]);
hold on;

% Add gray horizontal 0 line
plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on;
xline(10.2233, 'r--', 'LineWidth', 2);
% Darken axes
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% Set limits and ticks
xlim([0, 60]);  % ~0 to 5 seconds in minutes
ylim([-2, 7]);
xticks([0,30, 60]);  % 0s, 2.5s, 5s in minutes
%xticklabels({'0', '2.5', '5'});  % Optional: keep in seconds for clarity

yticks([-2,0, 7]);
yticklabels({'-2', '0', '7'});

% Aesthetics
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);

xlabel('time (minutes)');  % or 'Time (s)' if you keep ticks in seconds
ylabel('mean df/F (z-scored)');
legend off;
box off;
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below is pie charts of pos/neg/non cells between Ghrelin/Saline
% Below is pie charts of semaglutide and saline pos/neg/non cells
rng(6401) %100
% Settings
mice_ghrelin = {'A', 'B', 'F'};
mice_saline = {'A','D','F'};
sema_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin/Mouse_%s/%s_ghrelin.mat';
saline_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin/Mouse_%s/%s_saline.mat';
Fs = 15;  % sampling rate
nine_min = 9 * 60 * Fs;
n_shuffles = 1000;
post_start_frame = 10000;

% Initialize cell modulation counts
sema_counts = [0 0 0];   % [pos neg non]
saline_counts = [0 0 0]; % [pos neg non]

for i = 1:length(mice_ghrelin)
    mouse_ghrelin = mice_ghrelin{i};
    mouse_saline = mice_saline{i};

    % Load sema
    sema_file = sprintf(sema_path, mouse_ghrelin, mouse_ghrelin);
    load(sema_file, 'sumC');
    data = sumC;

    n_cells = size(data, 1);
    n_frames = size(data, 2);

    pre = mean(data(:, 1:nine_min), 2);
    post_real = mean(data(:, end-nine_min+1:end), 2);
    real_diff = post_real - pre;

    null_diffs = nan(n_cells, n_shuffles);
    for s = 1:n_shuffles
        idx = randi([post_start_frame, n_frames - nine_min + 1]);
        post_shuffled = mean(data(:, idx:idx+nine_min-1), 2);
        null_diffs(:, s) = post_shuffled - pre;
    end

    % Percentile thresholds
    lower_thresh = prctile(null_diffs, 15, 2);
    upper_thresh = prctile(null_diffs, 85, 2);

    is_neg = real_diff < lower_thresh;
    is_pos = real_diff > upper_thresh;
    is_non = ~is_neg & ~is_pos;

    sema_counts = sema_counts + [sum(is_pos), sum(is_neg), sum(is_non)];

    % Load saline
    saline_file = sprintf(saline_path, mouse_saline, mouse_saline);
    load(saline_file, 'sumC_sal');
    data = sumC_sal;

    n_cells = size(data, 1);
    n_frames = size(data, 2);

    pre = mean(data(:, 1:nine_min), 2);
    post_real = mean(data(:, end-nine_min+1:end), 2);
    real_diff = post_real - pre;

    null_diffs = nan(n_cells, n_shuffles);
    for s = 1:n_shuffles
        idx = randi([post_start_frame, n_frames - nine_min + 1]);
        post_shuffled = mean(data(:, idx:idx+nine_min-1), 2);
        null_diffs(:, s) = post_shuffled - pre;
    end

    lower_thresh = prctile(null_diffs, 15, 2);
    upper_thresh = prctile(null_diffs, 85, 2);

    is_neg = real_diff < lower_thresh;
    is_pos = real_diff > upper_thresh;
    is_non = ~is_neg & ~is_pos;

    saline_counts = saline_counts + [sum(is_pos), sum(is_neg), sum(is_non)];
end

% Plot pie charts
figure;

subplot(1,2,1);
pie(sema_counts, {'Positively Modulated', 'Negatively Modulated', 'Non-modulated'});
title('Ghrelin Group');

subplot(1,2,2);
pie(saline_counts, {'Positively Modulated', 'Negatively Modulated', 'Non-modulated'});%
title('Saline Group');
%exportgraphics(gca,'GhrelinPieCharts_Saline.eps','ContentType','vector')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Heatmap for Ghrelin
% === Set up ===
base_dir = '/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin';
mice = {'A', 'B', 'F'};
all_cells = [];

for i = 1:length(mice)
    mouse = mice{i};
    folder = fullfile(base_dir, ['Mouse_' mouse]);
    file = fullfile(folder, [mouse '_ghrelin.mat']);
    
    % Load sumC
    data = load(file, 'sumC'); % Change to sumC_sal for saline
    sumC = data.sumC; % Change to data.sumC_sal for saline
    
    % Check length and truncate
    if size(sumC, 2) >= 54000
        sumC = sumC(:, 1:54000);
    else
        warning('%s has fewer than 54000 frames, skipping.', mouse);
        continue
    end
    % Remove frames 8550 to 10350 (inclusive)
    sumC(:, 8550:10350) = [];
    
    % Concatenate
    all_cells = [all_cells; sumC];
end

% Sort by peak time in post-injection window (frame 9000 to 54000)

[~, avg_activity] = max(all_cells,[], 2);  % Average across time for each cell
[~, sort_idx] = sort(avg_activity, 'ascend');
sorted_data = all_cells(sort_idx, :);

% Plot
figure('Color','w');
imagesc(sorted_data);
%cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1);  % blue to white
%        ones(128,1), linspace(1,0,128)', linspace(1,0,128)']; % white to red
colormap('hot');
cb = colorbar;  % get handle to the colorbar
% Set ticks only at the min and max
clim([-1, 4]);
cb.Ticks = [-1, 4];

% Label as 'min' and 'max'
cb.TickLabels = {'min', 'max'};

% Optional: remove colorbar label (title)
cb.Label.String = '';
hold on;
xlabel('time (minutes)');
ylabel('cell #');
xlim([0 52200]);
yticks([1 size(sorted_data, 1)]);
xticks([0 26100 52200]);
xticklabels({'0','30','60'});
set(gca, 'TickDir','in','fontsize', 24)

% Add injection line
hold on;
x = 8550.5;  % injection frame at 10 min
xline(x, '--', 'ghrelin', ... % Change to Saline for saline heatmap
    'LineWidth', 3,...
    'Color', [1 1 1], ...       % white line
    'FontWeight', 'bold', ...
    'FontSize', 24, ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom');

exportgraphics(gca,'MS_GhrelinHeatmap_2minExcl.eps','ContentType','vector')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is example traces from saline and ghrelin injected mouse F
% Note: only 15 cells were used from the total 17 co-registered cells

load("/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin/Mouse_F/F_ghrelin.mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin/Mouse_F/F_saline(old).mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Ghrelin/Mouse_F/F_ghrelin_coreg(new2).mat");
cutoff = 1.5;    % Hz cutoff for Gaussian smoothing
frame_rate = 15;
% Gaussian filter kernel
sigma = frame_rate / (2 * pi * cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2 * sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

 % === Apply Gaussian filter ===
    filtered_sema = zeros(size(sumC));
    filtered_saline = zeros(size(sumC_sal));

    for j = 1:size(sumC, 1)
        filtered_leptin(j, :) = conv(sumC(j, :), gauss_kernel, 'same');
    end
    for j = 1:size(sumC_sal, 1)
        filtered_saline(j, :) = conv(sumC_sal(j, :), gauss_kernel, 'same');
    end




A = cell_registered_struct.cell_to_index_map; % CellReg output of co-registered indices
A(any(A == 0, 2), :) = [];


rows = A(:,1); % indices of leptin mice
rows_sal = A(:,2); % indices of corresponding saline mice

data = filtered_leptin(rows,:); % sumC is leptin traces from all cells from A_leptin.mat
data_sal = filtered_saline(rows_sal, :); % sumC_sal is saline traces from all cells from A_saline.mat
offset_factor = 5;
fps = 15;


data(7,:) = [];
data(9,:) = [];
data_sal(7,:) = [];
data_sal(9,:) = [];

% Set time vectors
time_sal = (1:size(data_sal, 2)) / fps;
time_leptin = (1:size(data, 2)) / fps;

% Determine y-limits based on data range + offsets
yl_sal = [0, size(data_sal, 1) * offset_factor + max(data_sal(:))];
yl_leptin = [0, size(data, 1) * offset_factor + max(data(:))];
[nCells, nFrames] = size(filtered_saline);
% Create tighter layout
figure;
%tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% Saline plot
%nexttile;
hold on
% Generate custom colormap: dark yellow → orange
endColor = [173, 206, 234] / 255;   % turquoise
startColor   = [93, 131, 174] / 255;  

cmap = zeros(nCells, 3);
for i = 1:3
    cmap(:, i) = linspace(startColor(i), endColor(i), nCells)';
end
for i = 1:size(data_sal, 1)
    plot(time_sal(52420:end), data_sal(i, 52420:end) + i * offset_factor, 'Color', cmap(i, :), 'LineWidth', 2);
end

xline(inject_sal/fps,"--r", "LineWidth",2);
xlabel('Time (seconds)');
ylabel('Amplitude');
axis off;
ylim(yl_sal);
xlim([52420/fps time_sal(end)]);
%exportgraphics(gca,'saline(ghrelin)_coreg_traces.eps','ContentType','vector')

% semaglutide plot

[nCells, nFrames] = size(filtered_sema);
%data = zscore(data, 0, 2);
% Create tighter layout
%tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% Saline plot
%nexttile;
hold on
% Generate custom colormap: dark yellow → orange
endColor = [148, 42, 80] / 255;  % dark yellow
startColor   = [198, 91, 118]  / 255; % red

cmap = zeros(nCells, 3);
for i = 1:3
    cmap(:, i) = linspace(startColor(i), endColor(i), nCells)';
end
figure;
hold on
for i = 1:size(data, 1)
    plot(time_leptin(54339:end), data(i, 54339:end) + i * offset_factor, 'Color', cmap(i, :), 'LineWidth', 2);
end
xline(9500/fps,"--r", "LineWidth",2);
axis off;

ylim(yl_leptin);
xlim([54339/fps time_leptin(end)]);
%exportgraphics(gca,'ghrelin_coreg_traces.eps','ContentType','vector')
