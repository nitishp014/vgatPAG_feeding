% Parameters
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding';
mice = {'Mouse_F'};
week = 'Week 2';
session = 'Session 2';
foods = {'Cheese', 'Chocolate', 'Broccoli', 'Carrot'};
food_labels_map = containers.Map({'Cheese', 'Chocolate', 'Broccoli', 'Carrot'}, 1:4);

for m = 1:length(mice)
    mouse = mice{m};
    fprintf('Processing %s...\n', mouse);

    mouse_letter = mouse(end);
    mat_file = fullfile(base_path, mouse, week, session, sprintf('%s_hungry.mat', mouse_letter));

    if ~isfile(mat_file)
        fprintf('  File not found: %s. Skipping %s.\n', mat_file, mouse);
        continue;
    end

    load(mat_file, 'sumC_hunger', 'timestamps');  % nCells x nFrames

    pca_data = [];
    food_labels = [];

    for f = 1:length(foods)
        food = foods{f};

        if ~isfield(timestamps.(mouse).Week2.Session2, food) || ...
           ~isfield(timestamps.(mouse).Week2.Session2.(food), 'Consumption')
            fprintf('  Skipping missing %s %s\n', mouse, food);
            continue;
        end

        vec = timestamps.(mouse).Week2.Session2.(food).Consumption;
        if isempty(vec), continue; end

        nBouts = floor(length(vec)/2);
        bouts = reshape(vec(1:2*nBouts), 2, [])';

        for b = 1:size(bouts,1)
            startF = bouts(b,1);
            endF = bouts(b,2);

            if isnan(startF) || isnan(endF) || endF > size(sumC_hunger,2)
                continue;
            end

            % Extract all frames for this bout
            bout_data = sumC_hunger(:, startF:endF)'; % frames x cells

            % Append to data and label each frame with the food label
            pca_data = [pca_data; bout_data];        % Each row = one frame's activity
            food_labels = [food_labels; repmat(food_labels_map(food), size(bout_data,1), 1)];
        end
    end

    if isempty(pca_data)
        fprintf('  No data for %s. Skipping PCA.\n', mouse);
        continue;
    end

    % Run PCA (choose 3 components for visualization)
    [coeff, score, ~, ~, explained] = pca(pca_data, 'NumComponents', 3);

    % Plot 3D PCA colored by food
    figure('Name', sprintf('3D PCA %s', mouse), 'NumberTitle', 'off');
    % Custom food-specific colors
colors = [
    227, 178, 69;  % Cheese - Dark Yellow
     103, 61, 37;  % Chocolate - Blue
    158, 209, 123;  % Broccoli - Green
    243, 163, 97   % Carrot - Purple
] / 255;  % Normalize to [0–1]
    hold on;
    for f = 1:length(foods)
        idx = food_labels == f;
        scatter3(score(idx,1), score(idx,2), score(idx,3), 10, colors(f,:), 'filled');
    end
    hold off;
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 30);
    xlabel('PC 1', 'FontSize', 30);
    ylabel('PC 2', 'FontSize', 30);
    zlabel('PC 3', 'FontSize', 30);
     % Determine axis limits
    x_min = min(score(:,1)); x_max = max(score(:,1));
    y_min = min(score(:,2)); y_max = max(score(:,2));
    z_min = min(score(:,3)); z_max = max(score(:,3));

    x_mid = (x_min + x_max) / 2;
    y_mid = (y_min + y_max) / 2;
    z_mid = (z_min + z_max) / 2;

    % Set ticks at min, mid, max
    xticks([x_min, x_mid, x_max]);
    yticks([y_min, y_mid, y_max]);
    zticks([z_min, z_mid, z_max]);

    % Set tick labels only for min and max
    xticklabels({num2str(x_min, '%.1f'), '', num2str(x_max, '%.1f')});
    yticklabels({num2str(y_min, '%.1f'), '', num2str(y_max, '%.1f')});
    zticklabels({num2str(z_min, '%.1f'), '', num2str(z_max, '%.1f')});
   % title(mouse);
    legend(foods, 'Location', 'northwestoutside');
    legend boxoff;
    view(45,30);
end
print(gcf, 'PCA_Foods.eps', '-depsc', '-painters');

% Parameters
base_path = '/Users/nitishpatel/Desktop/Miniscope Feeding';
mice = {'Mouse_A', 'Mouse_B', 'Mouse_D', 'Mouse_E', 'Mouse_F'};
week = 'Week 2';
session = 'Session 2';
foods = {'Cheese', 'Chocolate', 'Broccoli', 'Carrot'};
food_labels_map = containers.Map({'Cheese', 'Chocolate', 'Broccoli', 'Carrot'}, 1:4);

avg_silhouette_scores = nan(length(mice),1);

for m = 1:length(mice)
    mouse = mice{m};
    fprintf('Processing %s...\n', mouse);

    mouse_letter = mouse(end);
    mat_file = fullfile(base_path, mouse, week, session, sprintf('%s_hungry.mat', mouse_letter));

    if ~isfile(mat_file)
        fprintf('  File not found: %s. Skipping %s.\n', mat_file, mouse);
        continue;
    end

    load(mat_file, 'sumC_hunger');  % nCells x nFrames

    pca_data = [];
    food_labels = [];

    for f = 1:length(foods)
        food = foods{f};

        if ~isfield(timestamps.(mouse).Week2.Session2, food) || ...
           ~isfield(timestamps.(mouse).Week2.Session2.(food), 'Consumption')
            fprintf('  Skipping missing %s %s\n', mouse, food);
            continue;
        end

        vec = timestamps.(mouse).Week2.Session2.(food).Consumption;
        if isempty(vec), continue; end

        nBouts = floor(length(vec)/2);
        bouts = reshape(vec(1:2*nBouts), 2, [])';

        for b = 1:size(bouts,1)
            startF = bouts(b,1);
            endF = bouts(b,2);

            if isnan(startF) || isnan(endF) || endF > size(sumC_hunger,2)
                continue;
            end

            bout_data = sumC_hunger(:, startF:endF)'; % frames x cells

            pca_data = [pca_data; bout_data];
            food_labels = [food_labels; repmat(food_labels_map(food), size(bout_data,1), 1)];
        end
    end

    if isempty(pca_data)
        fprintf('  No data for %s. Skipping PCA and silhouette.\n', mouse);
        continue;
    end

    % Run PCA (3 components)
    [~, score, ~, ~, ~] = pca(pca_data, 'NumComponents', 3);

    % Compute silhouette values using the food labels as clusters
    sil_values = silhouette(score, food_labels);

    % Average silhouette score for this mouse
    avg_silhouette_scores(m) = mean(sil_values);

    fprintf('  %s average silhouette score: %.4f\n', mouse, avg_silhouette_scores(m));
end

% Remove NaNs before averaging
valid_scores = avg_silhouette_scores(~isnan(avg_silhouette_scores));
mean_sil_score = mean(valid_scores);

% Calculate mean and SEM
mean_sil = mean(valid_scores);
sem_sil = std(valid_scores) / sqrt(length(valid_scores));

% Plot
figure;
hold on;

% Bar for mean
%bar(0.15, mean_sil, 'FaceColor', [62/255, 45/255, 178/255], 'BarWidth', 0.15);
bar(0.15, mean_sil, 'FaceColor', [254/255, 216/255, 161/255], 'BarWidth', 0.15, 'EdgeColor','none');
% Error bar (SEM)
errorbar(0.15, mean_sil, sem_sil, 'k', 'LineWidth', 2, 'CapSize', 12);

% Scatter individual points
%scatter(repelem(0.15, length(valid_scores)), valid_scores, 50, [0.6,0.6 0.6], 'filled');
scatter(repelem(0.15, length(valid_scores)), valid_scores, 50, [0 0 0], 'filled');

% Y-axis formatting
ylim([-0.1 0.8]);
yticks([-0.1 0.8]);
yline(0,'r--', 'LineWidth', 2);

% X-axis
xlim([0 0.3]);


% Axis formatting
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 30, ...
         'Box', 'off', 'XColor', 'none');

ylabel({'clustering quality'; '(silhouette score)'}, 'FontSize', 30);

% Statistical test
[h, p] = ttest(valid_scores, 0); % two-tailed test
fprintf('One-sample t-test (silhouette > 0): h=%d, p=%.4f\n', h, p);


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

text(0.15, 0.85, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.15, 0.85, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

