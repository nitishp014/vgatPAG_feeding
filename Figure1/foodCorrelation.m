% Setup
mice = {'A','B','D','E','F','G'};
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding';
Fs = 15;
cutoff = 1.5;
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);
food_types = {'Cheese', 'Broccoli', 'Chocolate', 'Carrot'};

% Store per-cell mean activity per bout, per food across mice
% Structure: bout_means_cell.(food) = cell array with rows: each cell's [bout1_mean, bout2_mean]
bout_means_cell = struct();
for f = 1:length(food_types)
    bout_means_cell.(food_types{f}) = [];  
    % Each row: [bout1_mean, bout2_mean] for one cell from any mouse
end

% Also store bout1 means for each food for cross-food comparisons
% We'll keep these in a struct of cell arrays per mouse to ensure alignment of cells
bout1_means_per_mouse = struct();
for m = 1:length(mice)
    bout1_means_per_mouse.(mice{m}) = struct();
    for f = 1:length(food_types)
        bout1_means_per_mouse.(mice{m}).(food_types{f}) = [];
    end
end

% Process each mouse
for m = 1:length(mice)
    letter = mice{m};
    file_path = fullfile(base_path, ['Mouse_' letter], 'Week 2', 'Session 2', [letter '_hungry.mat']);
    if ~isfile(file_path)
        fprintf('Missing: %s\n', file_path);
        continue;
    end

    fprintf('Processing Mouse_%s...\n', letter);
    data = load(file_path);
    data_hunger = data.sumC_hunger;  % cells x time

    % Filter traces with Gaussian kernel
    filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end

    % For each food, extract bouts and calculate per-cell mean activity in each bout
    for f = 1:length(food_types)
        food = food_types{f};

        % Try to get consumption bouts timestamps for this food
        try
            bouts = getfield(data.timestamps, ['Mouse_' letter], 'Week2', 'Session2', food, 'Consumption');
        catch
            % No bouts for this food/mouse
            continue;
        end

        % Need at least 2 bouts (start,end pairs for 2 bouts)
        if length(bouts) < 4
            continue;
        end

        % Calculate per-cell mean activity in bouts
        % bouts is vector like [start1 end1 start2 end2 ...]
        bout1_indices = bouts(1):bouts(2);
        bout2_indices = bouts(3):bouts(4);

        bout1_mean_cells = mean(filtered(:, bout1_indices), 2); % mean activity per cell in bout1
        bout2_mean_cells = mean(filtered(:, bout2_indices), 2); % mean activity per cell in bout2

        % Append per-cell bout means to global struct for this food
        bout_means_cell.(food) = [bout_means_cell.(food); [bout1_mean_cells, bout2_mean_cells]];

        % Store bout1 means per cell for cross-food correlations, per mouse
        bout1_means_per_mouse.(letter).(food) = bout1_mean_cells;
    end
end

% Now calculate correlation matrix

n = length(food_types);
cross_corr = nan(n, n);
pval_matrix = nan(n, n);

% Diagonal: correlate bout1 vs bout2 mean across all cells pooled across mice for same food
for i = 1:n
    food = food_types{i};
    data_i = bout_means_cell.(food); % Nx2 matrix: bout1, bout2 per cell
    if size(data_i,1) >= 2
        [r, p] = corr(data_i(:,1), data_i(:,2), 'rows', 'complete');
        cross_corr(i,i) = r;
        pval_matrix(i,i) = p;
    end
end

% Off-diagonal: correlate bout1 means across cells between two different foods
% Must match cells by mouse and cell index
for i = 1:n
    for j = 1:n
        if i == j
            continue; % already done diagonal
        end

        food_i = food_types{i};
        food_j = food_types{j};

        % Collect vectors of matched cells' bout1 means across all mice
        vec_i = [];
        vec_j = [];

        for m = 1:length(mice)
            letter = mice{m};
            bi = bout1_means_per_mouse.(letter).(food_i);
            bj = bout1_means_per_mouse.(letter).(food_j);

            % Check if both exist and have same number of cells
            if isempty(bi) || isempty(bj)
                continue;
            end
            len = min(length(bi), length(bj));
            if len < 2
                continue; % too few cells for correlation
            end

            % Append matched cells' bout1 means for these two foods
            vec_i = [vec_i; bi(1:len)];
            vec_j = [vec_j; bj(1:len)];
        end

        if length(vec_i) >= 2
            [r, p] = corr(vec_i, vec_j, 'rows', 'complete');
            cross_corr(i,j) = r;
            pval_matrix(i,j) = p;
        end
    end
end

% Plotting
figure;
imagesc(cross_corr);
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {num2str(0), num2str(1)};
clim([0 1]);
axis equal tight;
xticks(1:n); 
yticks(1:n);
xticklabels(food_types); yticklabels(food_types);
set(gca, 'FontSize', 30);
%%%%%%%%%%%
% Bar graph for correlation matrix
% Extract diagonal and off-diagonal elements
diag_vals = diag(cross_corr);
off_diag_vals = cross_corr(tril(true(size(cross_corr)), -1));  % -1 excludes diagonal


% Calculate means
mean_diag = mean(diag_vals);
mean_off_diag = mean(off_diag_vals);

figure;
b = bar(0.15, [mean(diag_vals)], 'FaceColor', 'flat');
hold on;
c = bar(0.35, [mean(off_diag_vals)], 'FaceColor', 'flat');
b.CData = [174/255, 65/255, 121/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [249/255 218/255 120/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.35 0.15], [mean(off_diag_vals), mean(diag_vals)], [std(off_diag_vals)/sqrt(length(off_diag_vals)), std(diag_vals)/sqrt(length(diag_vals))], 'k.', 'LineWidth', 1);
hold on;
plot(0.35, off_diag_vals, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.15, diag_vals, 'ok', 'MarkerFaceColor', 'k');
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)

xlim([0 0.5])
yticks([0 1])
ylim([0 1])
legend({'diagonals', 'off-diagonals'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel('correlation (r)', 'FontSize', 24);
p = ranksum(diag_vals, off_diag_vals);
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

text(0.24, 1.05, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 1.08, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
exportgraphics(gca,'correlation_bars.eps','ContentType','vector')















% Above is correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Below is example correlations from Mouse B specifically (83 cells)!
% Collect all cells' bout1 and bout2 means for Cheese across mice
bout1_all_cells = [];
bout2_all_cells = [];

for m = 1:length(mice)
    letter = mice{m};
    file_path = fullfile(base_path, ['Mouse_' letter], 'Week 2', 'Session 2', [letter '_hungry.mat']);
    if ~isfile(file_path), continue; end

    data = load(file_path);
    data_hunger = data.sumC_hunger;

    % Filter traces with Gaussian kernel
    filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end

    % Get Cheese consumption bouts timestamps
    try
        bouts = getfield(data.timestamps, ['Mouse_' letter], 'Week2', 'Session2', 'Cheese', 'Consumption');
    catch
        continue;
    end

    if length(bouts) < 4, continue; end  % Need at least two bouts

    % Calculate mean activity per cell during each Cheese bout
    bout1_cell_means = mean(filtered(:, bouts(1):bouts(2)), 2); % cells x 1
    bout2_cell_means = mean(filtered(:, bouts(3):bouts(4)), 2);

    % Append to overall arrays
    bout1_all_cells = [bout1_all_cells; bout1_cell_means];
    bout2_all_cells = [bout2_all_cells; bout2_cell_means];
end

% Plot correlation scatter for Cheese bouts across all cells
figure;
scatter(bout1_all_cells, bout2_all_cells, 100, 'filled', 'MarkerFaceColor', [227/255 178/255 69/255]);
% Line of best fit
coeffs = polyfit(bout1_all_cells, bout2_all_cells, 1);  % linear fit
x_fit = linspace(-2, 3, 100);  % domain
y_fit = polyval(coeffs, x_fit);  % predicted values
hold on;
plot(x_fit, y_fit, '-', 'Color', [227/255 178/255 69/255], 'LineWidth', 3);  % same yellow as scatter
xlabel({'mean df/F (z-scored)'; 'cheese bout #1'});
ylabel({'mean df/F (z-scored)'; 'cheese bout #2'});
ylim([-2 3])
xlim([-2 3])
yticks([-2 3])
xticks([-2 3])
set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 30);
grid off;

% Calculate and display correlation coefficient on plot
[r, p] = corr(bout1_all_cells, bout2_all_cells, 'rows', 'complete');
text(-1, 2, sprintf('r = %.2f', r), 'HorizontalAlignment', 'center', 'FontSize', 30, 'FontWeight', 'bold');

% Initialize arrays to hold per-cell bout1 means for Cheese and Carrot across mice
cheese_bout1_all = [];
carrot_bout1_all = [];

for m = 1:length(mice)
    letter = mice{m};
    file_path = fullfile(base_path, ['Mouse_' letter], 'Week 2', 'Session 2', [letter '_hungry.mat']);
    if ~isfile(file_path), continue; end

    data = load(file_path);
    data_hunger = data.sumC_hunger;

    % Filter traces
    filtered = zeros(size(data_hunger));
    for j = 1:size(data_hunger, 1)
        filtered(j, :) = conv(data_hunger(j, :), gauss_kernel, 'same');
    end

    % Try to get Cheese bouts
    try
        cheese_bouts = getfield(data.timestamps, ['Mouse_' letter], 'Week2', 'Session2', 'Cheese', 'Consumption');
    catch
        cheese_bouts = [];
    end

    % Try to get Carrot bouts
    try
        carrot_bouts = getfield(data.timestamps, ['Mouse_' letter], 'Week2', 'Session2', 'Chocolate', 'Consumption');
    catch
        carrot_bouts = [];
    end

    % Need at least first bout for both foods
    if length(cheese_bouts) < 2 || length(carrot_bouts) < 2
        continue;
    end

    % Calculate mean activity per cell during Cheese bout 1 and Carrot bout 1
    cheese_bout1_cell_means = mean(filtered(:, cheese_bouts(1):cheese_bouts(2)), 2);
    carrot_bout1_cell_means = mean(filtered(:, carrot_bouts(1):carrot_bouts(2)), 2);

    % Only keep the minimum number of cells (if different number of cells, trim)
    n_cells = min(length(cheese_bout1_cell_means), length(carrot_bout1_cell_means));
    cheese_bout1_cell_means = cheese_bout1_cell_means(1:n_cells);
    carrot_bout1_cell_means = carrot_bout1_cell_means(1:n_cells);

    % Append
    cheese_bout1_all = [cheese_bout1_all; cheese_bout1_cell_means];
    carrot_bout1_all = [carrot_bout1_all; carrot_bout1_cell_means];
end

% Plot correlation scatter for Cheese bouts across all cells

   % 227/255 178/255 69/255;   Cheese - Dark Yellow
   %  103/255 61/255 37/255;   Chocolate - Brown
   % 158/255 209/255 123/255;   Broccoli - Green
   %  243/255 163/255 97/255    Carrot - Orange
figure;
scatter(cheese_bout1_all, carrot_bout1_all, 100, 'filled', 'MarkerFaceColor', [103/255 61/255 37/255]);
% Line of best fit
coeffs = polyfit(cheese_bout1_all, carrot_bout1_all, 1);  % linear fit
x_fit = linspace(-2, 3, 100);  % domain
y_fit = polyval(coeffs, x_fit);  % predicted values
hold on;
plot(x_fit, y_fit, '-', 'Color', [103/255 61/255 37/255], 'LineWidth', 3);  % same green as scatter

xlabel({'mean df/F (z-scored)'; 'cheese bout #1'});
ylabel({'mean df/F (z-scored)'; 'chocolate bout #1'});
ylim([-2 3])
xlim([-2 3])
yticks([-2 3])
xticks([-2 3])
set(gca,'linewidth',2, 'TickDir','out', 'fontsize', 30);
grid off;

% Calculate and display correlation coefficient on plot
[r, p] = corr(cheese_bout1_all, carrot_bout1_all, 'rows', 'complete');
text(2, 2, sprintf('r = %.2f', r), 'HorizontalAlignment', 'center', 'FontSize', 30, 'FontWeight', 'bold');
%change to 2, 2 for chocolate




