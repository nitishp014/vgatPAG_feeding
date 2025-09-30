% Define parameters
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Hunger_CoReg';
mouse_list = {'Mouse_A','Mouse_B','Mouse_D','Mouse_F','Mouse_G'}; %Note Mouse_E was no include in analysis
Fs = 15;               % Sampling frequency 
cutoff = 1.5;          % Cutoff frequency
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0
    win_size = win_size + 1;
end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

% Storage for all cell averages
all_sated_cell_means = [];
all_hunger_cell_means = [];

% Loop through mice
for i = 1:length(mouse_list)
    mouse_folder = fullfile(base_path, mouse_list{i});
    mat_files = dir(fullfile(mouse_folder, '*_CoReg.mat'));
    
    if isempty(mat_files)
        fprintf('No CoReg file found for %s\n', mouse_list{i});
        continue;
    end
    
    % Load the file
    load(fullfile(mouse_folder, mat_files(1).name), 'cell_registered_struct', 'sumC', 'sumC_hunger');
    
    A = cell_registered_struct.cell_to_index_map;
    A(any(A == 0, 2), :) = [];
    rows = A(:, 2); % sated
    rows_hungry = A(:, 1); % hungry

    % Extract and filter sated data
    data = sumC(rows, :);
    data_filtered = zeros(size(data));
    for j = 1:size(data, 1)
        data_filtered(j, :) = conv(data(j, :), gauss_kernel, 'same');
    end
    sated_means = mean(data_filtered(:, 1:1800), 2);  % per-cell mean
    all_sated_cell_means = [all_sated_cell_means; sated_means];

    % Extract and filter hunger data
    data_hunger = sumC_hunger(rows_hungry, :);
    data_hunger_filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        data_hunger_filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end
    hunger_means = mean(data_hunger_filtered(:, 1:1800), 2);  % per-cell mean
    all_hunger_cell_means = [all_hunger_cell_means; hunger_means];
end


% Compute means and SEMs
group_means = [mean(all_hunger_cell_means), mean(all_sated_cell_means)];
group_sems = [std(all_hunger_cell_means)/sqrt(length(all_hunger_cell_means)), ...
              std(all_sated_cell_means)/sqrt(length(all_sated_cell_means))];

% Bar plot
figure;
hold on;
c = bar(0.15, group_means(1), 'FaceColor', 'flat'); % hunger
b = bar(0.35, group_means(2), 'FaceColor', 'flat'); % sated

c.CData = [65/255 111/255 178/255];
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
b.CData = [249/255 218/255 120/255]; 
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
% Add error bars
errorbar([0.15 0.35], group_means, group_sems, 'k.', 'LineWidth', 1);

% Aesthetics
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
xlim([0 0.5])
yticks([-0.3 0.5])
ylim([-0.3 0.5]);
ylabel('mean df/F (z-scored)');
legend({'hunger', 'co-registered sated'}, 'Location', 'eastoutside');
legend boxoff
box off;

p = signrank(all_hunger_cell_means, all_sated_cell_means);

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

text(0.24, 0.543, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
text(0.24, 0.575, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

hold off;



% Above is for bar graph comparing hungry vs sated mice activity (2-min)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is for bar graph comparing hungry vs sated mice positive AUC
% (2-min)
% Define parameters
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding/Hunger_CoReg';
mouse_list = {'Mouse_A','Mouse_B','Mouse_D','Mouse_F','Mouse_G'};
Fs = 15;               % Sampling frequency
cutoff = 1.5;          % Cutoff frequency
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0
    win_size = win_size + 1;
end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);
timestep = 1 / Fs;

% Storage for all cell AUCs
all_sated_cell_aucs = [];
all_hunger_cell_aucs = [];

% Loop through mice
for i = 1:length(mouse_list)
    mouse_folder = fullfile(base_path, mouse_list{i});
    mat_files = dir(fullfile(mouse_folder, '*_CoReg.mat'));
    
    if isempty(mat_files)
        fprintf('No CoReg file found for %s\n', mouse_list{i});
        continue;
    end
    
    % Load the file
    load(fullfile(mouse_folder, mat_files(1).name), 'cell_registered_struct', 'sumC', 'sumC_hunger');
    
    A = cell_registered_struct.cell_to_index_map;
    A(any(A == 0, 2), :) = [];
    rows = A(:, 2); % sated
    rows_hungry = A(:, 1); % hungry

    % Extract and filter sated data
    data = sumC(rows, :);
    data_filtered = zeros(size(data));
    for j = 1:size(data, 1)
        data_filtered(j, :) = conv(data(j, :), gauss_kernel, 'same');
    end
    sated_aucs = zeros(size(data_filtered, 1), 1);
    for j = 1:size(data_filtered, 1)
        pos_vals = data_filtered(j, 1:1800);
        pos_vals = pos_vals(pos_vals > 0);
        sated_aucs(j) = sum(pos_vals) * timestep;
    end
    all_sated_cell_aucs = [all_sated_cell_aucs; sated_aucs];

    % Extract and filter hunger data
    data_hunger = sumC_hunger(rows_hungry, :);
    data_hunger_filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        data_hunger_filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end
    hunger_aucs = zeros(size(data_hunger_filtered, 1), 1);
    for j = 1:size(data_hunger_filtered, 1)
        pos_vals = data_hunger_filtered(j, 1:1800);
        pos_vals = pos_vals(pos_vals > 0);
        hunger_aucs(j) = sum(pos_vals) * timestep;
    end
    all_hunger_cell_aucs = [all_hunger_cell_aucs; hunger_aucs];
end

% Compute group means and SEMs
group_means = [mean(all_hunger_cell_aucs), mean(all_sated_cell_aucs)];
group_sems = [std(all_hunger_cell_aucs)/sqrt(length(all_hunger_cell_aucs)), ...
              std(all_sated_cell_aucs)/sqrt(length(all_sated_cell_aucs))];

% Bar plot
figure;
hold on;
c = bar(0.15, group_means(1), 'FaceColor', 'flat'); % hunger
b = bar(0.35, group_means(2), 'FaceColor', 'flat'); % sated

c.CData = [65/255 111/255 178/255];
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
b.CData = [249/255 218/255 120/255]; 
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
% Add error bars
errorbar([0.15 0.35], group_means, group_sems, 'k.', 'LineWidth', 1);

% Aesthetics
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
xlim([0 0.5])
yticks([0 100])
ylim([0 100]);
ylabel('positive AUC (z-scored df/F)');
legend({'hunger', 'co-registered sated'}, 'Location', 'eastoutside');
legend boxoff
box off;

% Significance test
p = signrank(all_hunger_cell_aucs, all_sated_cell_aucs);

% Asterisk notation
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

text(0.24, 104, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % significance
text(0.24, 107, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is example traces of hungry mice (from Mouse_A)
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_A/Week 2/Session 2/A_hungry.mat");
data = sumC_hunger(rows_hungry,1:1800);

% Smooth with Gaussian
Fs = 15;
cutoff = 2;
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

filtered = zeros(size(data));
for j = 1:size(data, 1)
    filtered(j, :) = conv(data(j, :), gauss_kernel, 'same');
end

% Offset for stacking
[nCells, nFrames] = size(filtered);
offset = 5;
offsets = (0:nCells-1)' * offset;
sumC_offset = filtered + offsets;


% Generate custom colormap: light blue → blue
endColor = [245, 189, 2] / 255;  % dark yellow
startColor   = [242, 131, 38]  / 255; % orange
cmap = zeros(nCells, 3);
for i = 1:3
    cmap(:, i) = linspace(startColor(i), endColor(i), nCells)';
end

% Plot each trace with its corresponding color
% Plot each trace with its corresponding color
figure; hold on;

timeVec = (0:nFrames-1) / 15 / 60;  % Convert frame indices to minutes

for i = 1:nCells
    plot(timeVec, sumC_offset(i, :), 'Color', cmap(i, :), 'LineWidth', 2);
end

xlim([-0.2, 2]);  % 0 to 2 minutes
ylim([-offset-0.2, offsets(end) + offset]);
xlabel('Time (min)');
ylabel('Fluorescence (offset)');

axis off;
% === Custom scale bars ===

% Constants
fps = 15;
xBarSeconds = 10;        % 30s horizontal scale bar
xBarMin = xBarSeconds / 60;
dfScale = 1;           % 40% dF/F in absolute trace units
cellIndex = 1;           % anchor to first trace
traceBase = offsets(cellIndex);  % vertical offset for that cell

% Choose origin (adjust xOrigin as needed)
xOrigin = 0;  % in minutes

% --- Horizontal time bar (30s) ---
line([xOrigin, xOrigin + xBarMin], [traceBase - 2.5, traceBase - 2.5], ...
     'Color', 'k', 'LineWidth', 3);
text(xOrigin + xBarMin/2, traceBase - 2.5, sprintf('%ds', xBarSeconds), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
%exportgraphics(gca,'MouseA_HungryTraces.eps','ContentType','vector')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is example traces of sated mice (from Mouse_A)
data = sumC(rows,1:1800);

% Smooth with Gaussian
Fs = 15;
cutoff = 2;
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

filtered = zeros(size(data));
for j = 1:size(data, 1)
    filtered(j, :) = conv(data(j, :), gauss_kernel, 'same');
end

% Offset for stacking
[nCells, nFrames] = size(filtered);
offset = 5;
offsets = (0:nCells-1)' * offset;
sumC_offset = filtered + offsets;

% Generate custom colormap: dark yellow → orange
endColor = [64, 224, 208] / 255;    % turquoise
startColor   = [79, 177, 178] / 255;   % slightly deeper light olive

cmap = zeros(nCells, 3);
for i = 1:3
    cmap(:, i) = linspace(startColor(i), endColor(i), nCells)';
end

% Plot each trace with its corresponding color
% Plot each trace with its corresponding color
figure; hold on;

timeVec = (0:nFrames-1) / 15 / 60;  % Convert frame indices to minutes

for i = 1:nCells
    plot(timeVec, sumC_offset(i, :), 'Color', cmap(i, :), 'LineWidth', 2);
end

xlim([-0.2, 2]);  % 0 to 2 minutes
ylim([-offset-0.2, offsets(end) + offset]);
xlabel('Time (min)');
ylabel('Fluorescence (offset)');

axis off;
% === Custom scale bars ===

% Constants
fps = 15;
xBarSeconds = 10;        % 30s horizontal scale bar
xBarMin = xBarSeconds / 60;
dfScale = 1;           % 40% dF/F in absolute trace units
cellIndex = 1;           % anchor to first trace
traceBase = offsets(cellIndex);  % vertical offset for that cell

% Choose origin (adjust xOrigin as needed)
xOrigin = 0;  % in minutes

% --- Horizontal time bar (30s) ---
line([xOrigin, xOrigin + xBarMin], [traceBase - 2.5, traceBase - 2.5], ...
     'Color', 'k', 'LineWidth', 3);
text(xOrigin + xBarMin/2, traceBase - 2.5, sprintf('%ds', xBarSeconds), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
exportgraphics(gca,'MouseA_SatedTraces.eps','ContentType','vector')
