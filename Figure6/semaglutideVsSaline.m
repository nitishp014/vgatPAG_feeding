% === Setup ===
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide';
mice = {'Mouse_A', 'Mouse_D', 'Mouse_F'};

Fs = 15;              % Sampling rate (Hz)
window_len = 8500;    % 9 min window (frames)
cutoff = 1.5;         % Gaussian filter cutoff (Hz)

% Gaussian filter
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% Storage
all_delta_sema = [];
all_delta_saline = [];

for i = 1:length(mice)
    mouse = mice{i};
    mouse_letter = mouse(end);
    mouse_path = fullfile(base_path, mouse);

    % File paths
    file_sema = fullfile(mouse_path, 'Semaglutide', sprintf('%s_sema.mat', mouse_letter));
    file_saline = fullfile(mouse_path, 'Saline', sprintf('%s_saline(old).mat', mouse_letter));
    file_coreg = fullfile(mouse_path, 'Semaglutide', sprintf('%s_sema_coreg(new3).mat', mouse_letter));

    % Skip missing files
    if ~isfile(file_sema) || ~isfile(file_saline) || ~isfile(file_coreg)
        fprintf('Skipping %s: missing file(s)\n', mouse);
        continue;
    end

    % Load data
    load(file_sema, 'sumC');                      % sema
    load(file_saline, 'sumC_sal', 'inject_sal');      % saline
    load(file_coreg, 'cell_registered_struct');   % CellReg

    A = cell_registered_struct.cell_to_index_map;
    A(any(A == 0, 2), :) = [];
    rows_sema = A(:,2);
    rows_saline = A(:,1);

    % Smooth
    filtered_sema = zeros(size(sumC));
    filtered_saline = zeros(size(sumC_sal));

    for j = 1:size(sumC,1)
        filtered_sema(j,:) = conv(sumC(j,:), gauss_kernel, 'same');
    end
    for j = 1:size(sumC_sal,1)
        filtered_saline(j,:) = conv(sumC_sal(j,:), gauss_kernel, 'same');
    end

    % Co-registered cells
    data_sema = filtered_sema(rows_sema, :);
    data_saline = filtered_saline(rows_saline, :);

    % --- SEMA windows ---
    inject_sema = 9450;
    pre_start_sema = 1;%floor((inject_sema - window_len)/2) + 1;
    pre_end_sema   = 10000;%pre_start_sema + window_len - 1;

    %post_length_sema = size(data_sema, 2) - inject_sema;
    post_start_sema  = 55000;%inject_sema + floor((post_length_sema - window_len)/2) + 1;
    post_end_sema    = 61000;%post_start_sema + window_len - 1;

    % --- SALINE windows ---
    inject_saline = inject_sal;
    pre_start_sal = 1900;%floor((inject_saline - window_len)/2) + 1;
    pre_end_sal   = 9800;%pre_start_sal + window_len - 1;

    post_length_sal = size(data_saline, 2) - inject_saline;
    post_start_sal  = 51000;%inject_saline + floor((post_length_sal - window_len)/2) + 1;
    post_end_sal    = 59000;%post_start_sal + window_len - 1;

    % --- Safety checks ---
    if any([pre_start_sema, pre_end_sema, post_start_sema, post_end_sema] > size(data_sema,2)) || ...
       any([pre_start_sal, pre_end_sal, post_start_sal, post_end_sal] > size(data_saline,2))
        fprintf('Skipping %s: window exceeds trace length.\n', mouse);
        continue;
    end

    % --- Compute deltas ---
    delta_sema = mean(data_sema(:, post_start_sema:end), 2) - ...
                 mean(data_sema(:, pre_start_sema:pre_end_sema), 2);

    delta_saline = mean(data_saline(:, post_start_sal:post_end_sal), 2) - ...
                   mean(data_saline(:, pre_start_sal:pre_end_sal), 2);

    all_delta_sema = [all_delta_sema; delta_sema];
    all_delta_saline = [all_delta_saline; delta_saline];
end

% === Signrank Test ===
if isempty(all_delta_sema) || isempty(all_delta_saline)
    error('No data available for statistical comparison.');
end

[p, ~, stats] = signrank(all_delta_sema, all_delta_saline);
fprintf('\nSignrank Test Results:\n');
fprintf('  p = %.4f\n', p);
fprintf('  z = %.2f\n', stats.zval);
fprintf('  N = %d co-registered cells\n', length(all_delta_sema));


figure;
b = bar(0.15, [mean(all_delta_sema)], 'FaceColor', 'flat');
hold on;
c = bar(0.35, [mean(all_delta_saline)], 'FaceColor', 'flat');
b.CData = [65/255 111/255 178/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [249/255 218/255 120/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.35 0.15], [mean(all_delta_saline), mean(all_delta_sema)], [std(all_delta_saline)/sqrt(length(all_delta_saline)), std(all_delta_sema)/sqrt(length(all_delta_sema))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 30)
%plot(0.15, ch_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
%plot(0.35, yfp_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
xlim([0 0.5])
yticks([-0.8 0.4])
ylim([-0.8 0.4])
legend({'semaglutide', 'saline'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'mean df/F (z-scored)'; '(post-pre)'}, 'FontSize', 30);


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

text(0.24, 0.45, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 0.49, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
%exportgraphics(gca,'MS_Updated_SemaglutideBars.eps','ContentType','vector')
% Heatmap below



% Heatmap code excluding injection-related activity (frames 8550 to 10350,
% a.k.a. mins (9:30 to 11:30)
% === Set up ===
base_dir = '/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide';
mice = {'A', 'D', 'F'};
all_cells = [];

for i = 1:length(mice)
    mouse = mice{i};
    folder = fullfile(base_dir, ['Mouse_' mouse]);
    folder = fullfile(folder, 'Semaglutide');
    file = fullfile(folder, [mouse '_sema.mat']);
    
    % Load sumC
    data = load(file, 'sumC');
    sumC = data.sumC;
    
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

% Sort by peak time in post-injection window (new frame indices!)
% Injection used to be at frame 9000, but now its between frames 8550 and
% 8551
injection_frame = 8550.5;  % 7199
[~, peak_activity] = max(all_cells,[], 2);
[~, sort_idx] = sort(peak_activity, 'ascend');  % or 'ascend' if preferred
sorted_data = all_cells(sort_idx, :);

% Plot
figure('Color','w');
imagesc(sorted_data);
%cmap = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1);  % blue to white
%        ones(128,1), linspace(1,0,128)', linspace(1,0,128)']; % white to red
colormap('hot');
cb = colorbar;  % get handle to the colorbar
clim([-1, 4]);
cb.Ticks = [-1, 4];
cb.TickLabels = {'min', 'max'};
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
xline(injection_frame, '--', 'Semaglutide', ...
    'LineWidth', 3,...
    'Color', [1 1 1], ...
    'FontWeight', 'bold', ...
    'FontSize', 24, ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom');

exportgraphics(gca,'MS_SemaglutideHeatmap_2minExcl.eps','ContentType','vector')




% Heatmap for Leptin (or Ghrelin, use A, B, F, G mice)
% === Set up ===
    base_dir = '/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin';
mice = {'A', 'D', 'F'};
all_cells = [];

for i = 1:length(mice)
    mouse = mice{i};
    folder = fullfile(base_dir, ['Mouse_' mouse]);
    file = fullfile(folder, [mouse '_saline.mat']);
    
    % Load sumC
    data = load(file, 'sumC_sal'); % Change to sumC_sal for saline
    sumC_sal = data.sumC_sal; % Change to data.sumC_sal for saline
    
    % Check length and truncate
    if size(sumC_sal, 2) >= 54000
        sumC_sal = sumC_sal(:, 1:54000);
    else
        warning('%s has fewer than 54000 frames, skipping.', mouse);
        continue
    end
    % Remove frames 8550 to 10350 (inclusive)
    sumC_sal(:, 8550:10350) = [];
    
    % Concatenate
    all_cells = [all_cells; sumC_sal];
end

% Sort by peak time in post-injection window (frame 9000 to 54000)

[~, peak_activity] = max(all_cells,[], 2);  % Average across time for each cell
[~, sort_idx] = sort(peak_activity, 'ascend');
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
xline(x, '--', 'Saline', ... % Change to Saline for saline heatmap
    'LineWidth', 3,...
    'Color', [1 1 1], ...       % white line
    'FontWeight', 'bold', ...
    'FontSize', 24, ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom');
%exportgraphics(gca,'MS_Saline(Sema)Heatmap_2minExcl.eps','ContentType','vector')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is correlation between co-registered cells from sema and saline
% (note we used the same saline that we used for Leptin)
% === Parameters ===
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide';
mice = {'Mouse_A','Mouse_D','Mouse_F'};
frame_rate = 15; % Hz
cutoff = 1.5;    % Hz cutoff for Gaussian smoothing
window_min = 9;
window_frames = frame_rate * 60 * window_min;

% Gaussian filter kernel
sigma = frame_rate / (2 * pi * cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2 * sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% Initialize containers
all_delta_sema = [];
all_delta_saline = [];

for i = 1:length(mice)
    mouse = mice{i};
    mouse_letter = mouse(end); % Assumes format "Mouse_X"
    mouse_path = fullfile(base_path, mouse);
    
    % Load .mat files
    sema_data = load(fullfile(mouse_path, 'Semaglutide', [mouse_letter '_sema.mat']));
    saline_data = load(fullfile(mouse_path, 'Saline', [mouse_letter '_saline(old).mat']));
    reg_data = load(fullfile(mouse_path, 'Semaglutide', [mouse_letter '_sema_coreg(new3).mat'])); % Assumes same name in each folder

    % Extract variables
    sumC = sema_data.sumC;
    inject = 9450;
    sumC_sal = saline_data.sumC_sal;
    inject_sal = saline_data.inject_sal;
    A = reg_data.cell_registered_struct.cell_to_index_map;

    % Remove rows with unregistered cells
    A(any(A == 0, 2), :) = [];
    rows_sema = A(:, 2); % indices into sumC
    rows_saline = A(:, 1); % indices into sumC_sal

    % === Apply Gaussian filter ===
    filtered_sema = zeros(size(sumC));
    filtered_saline = zeros(size(sumC_sal));

    for j = 1:size(sumC, 1)
        filtered_sema(j, :) = conv(sumC(j, :), gauss_kernel, 'same');
    end
    for j = 1:size(sumC_sal, 1)
        filtered_saline(j, :) = conv(sumC_sal(j, :), gauss_kernel, 'same');
    end

    % Only keep co-registered cells
    data_sema = filtered_sema(rows_sema, :);
    data_saline = filtered_saline(rows_saline, :);

    % === Pre- and post-window indexing ===
    pre_idx_sema = max(1, inject - window_frames) : inject - 1;
    post_idx_sema = 52000:60000; % For correlation plot: use leptin 52000:60000
                                   % For bar graph: use leptin 20000:30000
    pre_idx_saline = max(1, inject_sal - window_frames) : inject_sal - 1;
    post_idx_saline = 52000:60000; % For correlation/bar plot: use saline 39000:60000

    % === Compute post - pre activity ===
    pre_sema = mean(data_sema(:, pre_idx_sema), 2);
    post_sema = mean(data_sema(:, post_idx_sema), 2);
    delta_sema = post_sema - pre_sema;

    pre_saline = mean(data_saline(:, pre_idx_saline), 2);
    post_saline = mean(data_saline(:, post_idx_saline), 2);
    delta_saline = post_saline - pre_saline;

    % === Accumulate results ===
    all_delta_sema = [all_delta_sema; delta_sema];
    all_delta_saline = [all_delta_saline; delta_saline];
end

% === Plot correlation ===
figure;
scatter(all_delta_sema, all_delta_saline, 90, 'filled', "MarkerFaceColor", [174/255, 65/255, 121/255],"LineWidth",1.2, "MarkerEdgeAlpha", 0.8);
xlabel({'semaglutide df/F (z-scored)'; '(post-pre)'});
ylabel({'saline df/F (z-scored)'; '(post-pre)'});

pfit = polyfit(all_delta_sema, all_delta_saline, 1);  % Returns [slope, intercept]
% Define x values using the axis limits
x_vals = xlim;  % or [min(ghrelin_mod), max(ghrelin_mod)] if you prefer
y_vals = polyval(pfit, x_vals);  % Evaluate the line at these x values
hold on;
% Plot the custom regression line
plot(x_vals, y_vals, '-', 'Color', [174/255, 65/255, 121/255], 'LineWidth', 1.5);

set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 24)
xlim([-3 3])
xticks([-3 0 3])
yticks([-2.5 0 2.5])
ylim([-2.5 2.5])

grid off;
box off;
%axis equal;
%lsline; % least-squares fit line

% === Correlation coefficient ===
[r, p] = corr(all_delta_sema, all_delta_saline);
text(0.5, 2.1, sprintf('r = %.2f\np = %.2f', r, p), ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'center');
%exportgraphics(gca,'MS_sema_saline_correlation.eps','ContentType','vector')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is example traces from saline and semaglutide injected mouse D
% Note: only 15 cells were used from the total 25 co-registered cells

load("/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_D/Semaglutide/D_sema.mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_D/Saline/D_saline(old).mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_D/Semaglutide/D_sema_coreg(new3).mat");
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


rows = A(:,2); % indices of leptin mice
rows_sal = A(:,1); % indices of corresponding saline mice

data = filtered_leptin(rows,:); % sumC is leptin traces from all cells from A_leptin.mat
data_sal = filtered_saline(rows_sal, :); % sumC_sal is saline traces from all cells from A_saline.mat
offset_factor = 5;
fps = 15;


data(18:24,:) = [];
data(13:15,:) = [];
data_sal(18:24,:) = [];
data_sal(13:15,:) = [];

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
%exportgraphics(gca,'saline(sema)_coreg_traces.eps','ContentType','vector')

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
%exportgraphics(gca,'sema_coreg_traces.eps','ContentType','vector')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is co-registration footprints from Mouse A
clear;
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_A/Semaglutide/A_sema_coreg(new3).mat");
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin/Mouse_A/Saline/footprints.mat");
A1 = permute(A, [2, 3, 1]);  % Now A is [280, 300, 90] Saline

load("/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_A/Semaglutide/footprints.mat");
A2 = permute(A, [2, 3, 1]);  % Now A is [250, 310, 61] Semaglutide

% === Parameters ===
sigma = 2;
threshold = 0.2;
alpha_val = 0.6;
[H, W, ~] = size(A1);

% === Colors ===
yellow = [249, 218, 120]/255;
blue   = [65, 111, 178]/255;
white  = [0.3, 0.3, 0.3];

% === CellReg Mapping ===
map = cell_registered_struct.cell_to_index_map;
co_registered_idx = find(map(:,1) > 0 & map(:,2) > 0);  % co-registered only
n_coreg = numel(co_registered_idx);
% === Pastel spectrum: red → violet in soft tones ===
% === Pastel spectrum: red → violet ===
hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb

% === Helper function ===
normalize = @(x) x / max(x(:));
cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);


% === ----- Plot A1 Figure ----- === THIS IS MOUSE D, WEEK 2, SESSION 2
overlay_rgb_A1 = zeros(H, W, 3);
overlay_alpha_A1 = zeros(H, W);

coreg_counter = 1;  % Index for cmap

for i = 1:size(map,1)
    idx1 = map(i, 1);
    idx2 = map(i, 2);

    if idx1 == 0, continue; end

    img = normalize(imgaussfilt(A1(:,:,idx1), sigma));
    mask = img > threshold;

    if idx2 > 0
        color = cmap(coreg_counter, :);  % Use next color from cmap
        coreg_counter = coreg_counter + 1;
    else
        color = white;  % A1-only
    end

    for ch = 1:3
        channel = overlay_rgb_A1(:,:,ch);
        channel(mask) = channel(mask) + color(ch);
        overlay_rgb_A1(:,:,ch) = channel;
    end
    overlay_alpha_A1(mask) = overlay_alpha_A1(mask) + alpha_val;
end

overlay_alpha_A1 = min(overlay_alpha_A1, 1);

figure('Color', 'w');
imshow(overlay_rgb_A1);
hold on;
h = imshow(overlay_rgb_A1);
set(h, 'AlphaData', overlay_alpha_A1);
axis image off;
%exportgraphics(gca,'MouseA_Saline(Sema)Footprints.eps','ContentType','vector')


% === Parameters ===
sigma = 2;
threshold = 0.2;
alpha_val = 0.6;
[H, W, ~] = size(A2);

% === Colors ===
yellow = [249, 218, 120]/255;
blue   = [65, 111, 178]/255;
white  = [0.3, 0.3, 0.3];

% === CellReg Mapping ===
map = cell_registered_struct.cell_to_index_map;
co_registered_idx = find(map(:,1) > 0 & map(:,2) > 0);  % co-registered only
n_coreg = numel(co_registered_idx);
% === Pastel spectrum: red → violet in soft tones ===
% === Pastel spectrum: red → violet ===
hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb

% Convert HSV to pastel RGB
cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);

% === Helper function ===
normalize = @(x) x / max(x(:));

% === ----- Plot A2 Figure ----- === THIS IS MOUSE D, WEEK 1, SESSION 3
overlay_rgb_A2 = zeros(H, W, 3);
overlay_alpha_A2 = zeros(H, W);

coreg_counter = 1;  % Index for cmap

for i = 1:size(map,1)
    idx1 = map(i, 1);
    idx2 = map(i, 2);

    if idx2 == 0, continue; end  % skip if not present in A2

    img = normalize(imgaussfilt(A2(:,:,idx2), sigma));
    mask = img > threshold;

    if idx1 > 0
        color = cmap(coreg_counter, :);  % Use next color from cmap
        coreg_counter = coreg_counter + 1;
    else
        color = white ;  % A2-only
    end

    for ch = 1:3
        channel = overlay_rgb_A2(:,:,ch);
        channel(mask) = channel(mask) + color(ch);
        overlay_rgb_A2(:,:,ch) = channel;
    end
    overlay_alpha_A2(mask) = overlay_alpha_A2(mask) + alpha_val;
end

overlay_alpha_A2 = min(overlay_alpha_A2, 1);

figure('Color', 'w');
imshow(overlay_rgb_A2);
hold on;
h = imshow(overlay_rgb_A2);
set(h, 'AlphaData', overlay_alpha_A2);
axis image off;
%exportgraphics(gca,'MouseA_SemaglutideFootprints.eps','ContentType','vector')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is pie charts of semaglutide and saline pos/neg/non cells
rng(6401) %100
% Settings
mice = {'A', 'D', 'F'};
sema_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Semaglutide/Mouse_%s/Semaglutide/%s_sema.mat';
saline_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin/Mouse_%s/%s_saline.mat';
Fs = 15;  % sampling rate
nine_min = 9 * 60 * Fs;
n_shuffles = 1000;
post_start_frame = 10000;

% Initialize cell modulation counts
sema_counts = [0 0 0];   % [pos neg non]
saline_counts = [0 0 0]; % [pos neg non]

for i = 1:length(mice)
    mouse = mice{i};

    % Load sema
    sema_file = sprintf(sema_path, mouse, mouse);
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
    saline_file = sprintf(saline_path, mouse, mouse);
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
title('Semaglutide Group');

subplot(1,2,2);
pie(saline_counts, {'Positively Modulated', 'Negatively Modulated', 'Non-modulated'});
%title('Saline Group');
%exportgraphics(gca,'SemaPieCharts_Sema.eps','ContentType','vector')

