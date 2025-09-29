A = cell_registered_struct.cell_to_index_map; % CellReg output of co-registered indices
A(any(A == 0, 2), :) = [];


rows = A(:,1); % indices of leptin mice
rows_sal = A(:,2); % indices of corresponding saline mice

data = filtered_leptin(rows,:); % sumC is leptin traces from all cells from A_leptin.mat
data_sal = filtered_saline(rows_sal, :); % sumC_sal is saline traces from all cells from A_saline.mat
offset_factor = 5;
fps = 15;

% Set time vectors
time_sal = (1:size(data_sal, 2)) / fps;
time_leptin = (1:size(data, 2)) / fps;

% Determine y-limits based on data range + offsets
yl_sal = [0, size(data_sal, 1) * offset_factor + max(data_sal(:))];
yl_leptin = [0, size(data, 1) * offset_factor + max(data(:))];

% Create tighter layout
figure;
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% leptin plot
nexttile;
hold on
for i = 1:size(data, 1)
    plot(time_leptin, data(i, :) + i * offset_factor, 'k', 'LineWidth', 1);
end
xline(inject/fps,"--r", "LineWidth",2);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Leptin');
ylim(yl_leptin);
xlim([0 time_leptin(end)]);

% Saline plot
nexttile;
hold on
for i = 1:size(data_sal, 1)
    plot(time_sal, data_sal(i, :) + i * offset_factor, 'k', 'LineWidth', 1);
end
xline(inject_sal/fps,"--r", "LineWidth",2);
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Saline');
ylim(yl_sal);
xlim([0 time_sal(end)]);
% Above is for plotting the traces of co-registered cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define paths
%To run this code: you need [Letter]_leptin_coreg.mat file with
%inject and inject_sal as the frame number at
%which injection occurred (can find this using the code above)

% === Parameters ===
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin';
mice = {'Mouse_D','Mouse_F','Mouse_G'};
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
all_delta_leptin = [];
all_delta_saline = [];

for i = 1:length(mice)
    mouse = mice{i};
    mouse_letter = mouse(end); % Assumes format "Mouse_X"
    mouse_path = fullfile(base_path, mouse);
    
    % Load .mat files
    leptin_data = load(fullfile(mouse_path, [mouse_letter '_leptin.mat']));
    saline_data = load(fullfile(mouse_path, [mouse_letter '_saline.mat']));
    reg_data = load(fullfile(mouse_path, [mouse_letter '_leptin_coreg.mat'])); % Assumes same name in each folder

    % Extract variables
    sumC = leptin_data.sumC;
    inject = leptin_data.inject;
    sumC_sal = saline_data.sumC_sal;
    inject_sal = saline_data.inject_sal;
    A = reg_data.cell_registered_struct.cell_to_index_map;

    % Remove rows with unregistered cells
    A(any(A == 0, 2), :) = [];
    rows_leptin = A(:, 1); % indices into sumC
    rows_saline = A(:, 2); % indices into sumC_sal

    % === Apply Gaussian filter ===
    filtered_leptin = zeros(size(sumC));
    filtered_saline = zeros(size(sumC_sal));

    for j = 1:size(sumC, 1)
        filtered_leptin(j, :) = conv(sumC(j, :), gauss_kernel, 'same');
    end
    for j = 1:size(sumC_sal, 1)
        filtered_saline(j, :) = conv(sumC_sal(j, :), gauss_kernel, 'same');
    end

    % Only keep co-registered cells
    data_leptin = filtered_leptin(rows_leptin, :);
    data_saline = filtered_saline(rows_saline, :);

    % === Pre- and post-window indexing ===
    pre_idx_leptin = max(1, inject - window_frames) : inject - 1;
    post_idx_leptin = 20000:30000; % For correlation plot: use leptin 52000:60000
                                   % For bar graph: use leptin 20000:30000
    pre_idx_saline = max(1, inject_sal - window_frames) : inject_sal - 1;
    post_idx_saline = 39000:60000; % For correlation/bar plot: use saline 39000:60000

    % === Compute post - pre activity ===
    pre_leptin = mean(data_leptin(:, pre_idx_leptin), 2);
    post_leptin = mean(data_leptin(:, post_idx_leptin), 2);
    delta_leptin = post_leptin - pre_leptin;

    pre_saline = mean(data_saline(:, pre_idx_saline), 2);
    post_saline = mean(data_saline(:, post_idx_saline), 2);
    delta_saline = post_saline - pre_saline;

    % === Accumulate results ===
    all_delta_leptin = [all_delta_leptin; delta_leptin];
    all_delta_saline = [all_delta_saline; delta_saline];
end

% === Plot correlation ===
figure;
scatter(all_delta_leptin, all_delta_saline, 90, 'filled', "MarkerFaceColor", [174/255, 65/255, 121/255],"LineWidth",1.2, "MarkerEdgeAlpha", 0.8);
xlabel({'leptin df/F (z-scored)'; '(post-pre)'});
ylabel({'saline df/F (z-scored)'; '(post-pre)'});

pfit = polyfit(all_delta_leptin, all_delta_saline, 1);  % Returns [slope, intercept]
% Define x values using the axis limits
x_vals = xlim;  % or [min(ghrelin_mod), max(ghrelin_mod)] if you prefer
y_vals = polyval(pfit, x_vals);  % Evaluate the line at these x values
hold on;
% Plot the custom regression line
plot(x_vals, y_vals, '-', 'Color', [174/255, 65/255, 121/255], 'LineWidth', 1.5);

set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 24)
xlim([-2.5 1.5])
xticks([-2.5 0 1.5])
yticks([-2.5 0 1.5])
ylim([-2.5 1.5])

grid off;
box off;
%axis equal;
%lsline; % least-squares fit line

% === Correlation coefficient ===
[r, p] = corr(all_delta_leptin, all_delta_saline);
text(0.5, -1.1, sprintf('r = %.2f\np = %.3f', r, p), ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'center');

%exportgraphics(gca,'leptin_saline_correlation.eps','ContentType','vector')



figure;
b = bar(0.15, [mean(all_delta_leptin)], 'FaceColor', 'flat');
hold on;
c = bar(0.35, [mean(all_delta_saline)], 'FaceColor', 'flat');
b.CData = [174/255, 65/255, 121/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [249/255 218/255 120/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.35 0.15], [mean(all_delta_saline), mean(all_delta_leptin)], [std(all_delta_saline)/sqrt(length(all_delta_saline)), std(all_delta_leptin)/sqrt(length(all_delta_leptin))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
%plot(0.15, ch_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
%plot(0.35, yfp_on_off, 'ok', 'MarkerFaceColor', 'k');
%hold on;
xlim([0 0.5])
yticks([-0.8 0.1])
ylim([-0.8 0.1])
legend({'leptin', 'saline'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'mean df/F (z-scored)'; '(post-pre)'}, 'FontSize', 24);
%exportgraphics(gca,'MS_leptin_saline_bars.eps','ContentType','vector')

signrank(all_delta_leptin, all_delta_saline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Below is for plotting Leptin #27 co-registered cell and Saline cell #12 (from sumC_sal) from Mouse_F as example traces
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin/Mouse_F/F_leptin.mat");
figure;
subplot(2, 1, 1);
avg_activity = sumC(27,:);
avg_activity = avg_activity(:, 1:54000);
frame_rate = 15;  % frames per second

% Convert x-values (frame indices) to time in minutes
time_minutes = (1:80:length(avg_activity)) / frame_rate / 60;

% Plot
plot(time_minutes, avg_activity(1:80:end), '-', 'LineWidth', 2, 'Color',[174/255, 65/255, 121/255]);
hold on;

% Add gray horizontal 0 line
plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on;
xline(10.1344, 'r--', 'LineWidth', 2);

% Darken axes
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% Set limits and ticks
xlim([0, 60]);  % ~0 to 5 seconds in minutes
ylim([-2, 7.5]);
xticks([0,30, 60]);  % 0s, 2.5s, 5s in minutes
%xticklabels({'0', '2.5', '5'});  % Optional: keep in seconds for clarity

yticks([-2,0, 7.5]);
yticklabels({'-2', '0', '7.5'});

% Aesthetics
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);

xlabel('time (minutes)');  % or 'Time (s)' if you keep ticks in seconds
ylabel('mean df/F (z-scored)');
legend off;
box off;
hold off;

load("/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin/Mouse_F/F_saline.mat");
data_sal = sumC_sal;
subplot(2, 1, 2);
avg_saline = data_sal(12,:);
time_saline = (1:80:length(avg_saline)) / frame_rate / 60;
plot(time_saline, avg_saline(1:80:end), '-', 'LineWidth', 2, 'Color',[231/255 173/255 82/255]);
hold on;

% Add gray horizontal 0 line
plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);
hold on;
xline(10.7567, 'r--', 'LineWidth', 2);
% Darken axes
ax = gca;
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;

% Set limits and ticks
xlim([0, 60]);  % ~0 to 5 seconds in minutes
ylim([-2, 7.5]);
xticks([0,30, 60]);  % 0s, 2.5s, 5s in minutes
%xticklabels({'0', '2.5', '5'});  % Optional: keep in seconds for clarity

yticks([-2,0, 7.5]);
yticklabels({'-2', '0', '7.5'});

% Aesthetics
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);

xlabel('time (minutes)');  % or 'Time (s)' if you keep ticks in seconds
ylabel('mean df/F (z-scored)');
legend off;
box off;
hold off;

%%%%%%%%%%%%%%%%%%%%%
% Code for Leptin heatmap
base_dir = '/Users/nitishpatel/Desktop/Miniscope Feeding/Leptin';
mice = {'A','D', 'F','G'};
all_cells = [];

for i = 1:length(mice)
    mouse = mice{i};
    folder = fullfile(base_dir, ['Mouse_' mouse]);
    %folder = fullfile(folder, 'Leptin');
    file = fullfile(folder, [mouse '_leptin.mat']);
    
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
xlabel('time from injection (min)');
ylabel('neuron #');
xlim([0 52200]);
yticks([1 size(sorted_data, 1)]);
xticks([0 26100 52200]);
xticklabels({'-10','20','50'});
set(gca, 'TickDir','in','fontsize', 24)

% Add injection line
xline(injection_frame, '--', 'leptin', ...
    'LineWidth', 3,...
    'Color', [1 1 1], ...
    'FontWeight', 'bold', ...
    'FontSize', 24, ...
    'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom');

%exportgraphics(gca,'MS_LeptinHeatmap_2minExcl.eps','ContentType','vector')