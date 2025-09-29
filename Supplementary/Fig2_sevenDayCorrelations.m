%% Parameters
baseDir = '/Users/nitishpatel/Desktop/Miniscope Feeding/Seven Days';
mice = {'A','D','F'};
days = [1,2,3,4,5,6,7];
Fs = 15;              % Hz
cutoff = 1.5;         % Gaussian filter cutoff
sigma = Fs / (2*pi*cutoff);
win_size = ceil(sigma * 6);
if mod(win_size, 2) == 0, win_size = win_size + 1; end
x = linspace(-win_size/2, win_size/2, win_size);
gauss_kernel = exp(-x.^2 / (2*sigma^2));
gauss_kernel = gauss_kernel / sum(gauss_kernel);

%% Load timestamps (eating bouts)
load(fullfile(baseDir, 'timestamps.mat')); % gives timestamps.Mouse_X.DayY

%% Prepare storage for pooled data
comparisons = {[1, 2],[1 3], [1 4],[1 5],[1 6],[1 7],[2 3],[2 4],[2 5],[2 6],[2 7],[3 4],[3 5],[3 6],[3 7],[4 5],[4 6],[4 7],[5 6],[5,7],[6 7]};
pooled = cell(3,2);  % {comparisonIndex, {x,y}}

for c = 1:numel(comparisons)
    pooled{c,1} = []; % x-values
    pooled{c,2} = []; % y-values
end

%% Helper function
toFrames = @(minSec,Fs) round(minSec(:,1)*60*Fs + minSec(:,2)*Fs + 1);

%% Loop through mice
for m = 1:numel(mice)
    mouse = mice{m};
    fprintf('Processing Mouse %s...\n', mouse);

    % Load and smooth sumC for each day
    sumC = struct();
    for d = 1:numel(days)
        day = days(d);
        dataFile = fullfile(baseDir, ...
            sprintf('Mouse_%s/Day %d/%s_day%d.mat', mouse, day, mouse, day));
        tmp = load(dataFile, 'sumC');

        % Smooth traces
        filtered = zeros(size(tmp.sumC));
        for j = 1:size(tmp.sumC,1)
            filtered(j,:) = conv(tmp.sumC(j,:), gauss_kernel, 'same');
        end
        sumC.(sprintf('day%d',day)) = filtered;
    end

    % Do each comparison
    for c = 1:numel(comparisons)
        d1 = comparisons{c}(1);
        d2 = comparisons{c}(2);

        % Load coregistration file
        regFile = fullfile(baseDir, sprintf('%s_%d-%d.mat', mouse, d1, d2));
        regStruct = load(regFile);
        A = regStruct.cell_registered_struct.cell_to_index_map;
        A(any(A==0,2),:) = [];  % remove unregistered

        rows_day1 = A(:,1); 
        rows_day2 = A(:,2);

        % --- Get pistachio bouts ---
        bouts_day1 = timestamps.(sprintf('Mouse_%s',mouse)).(sprintf('Day%d',d1));
        bouts_day2 = timestamps.(sprintf('Mouse_%s',mouse)).(sprintf('Day%d',d2));

        % Convert bouts to frames
        startF1 = toFrames(bouts_day1(:,1:2), Fs);
        endF1   = toFrames(bouts_day1(:,3:4), Fs);
        startF2 = toFrames(bouts_day2(:,1:2), Fs);
        endF2   = toFrames(bouts_day2(:,3:4), Fs);

        % Clamp to recording length
        nFrames1 = size(sumC.(sprintf('day%d',d1)),2);
        nFrames2 = size(sumC.(sprintf('day%d',d2)),2);

        startF1 = max(1, min(startF1, nFrames1));
        endF1   = max(1, min(endF1,   nFrames1));
        startF2 = max(1, min(startF2, nFrames2));
        endF2   = max(1, min(endF2,   nFrames2));

        % --- Compute averages ---
        avg_day1 = nan(length(rows_day1),1);
        avg_day2 = nan(length(rows_day2),1);

        for k = 1:length(rows_day1)
            trace1 = sumC.(sprintf('day%d',d1))(rows_day1(k),:);
            trace2 = sumC.(sprintf('day%d',d2))(rows_day2(k),:);

            vals1 = [];
            for b = 1:numel(startF1)
                if endF1(b) > startF1(b)
                    vals1 = [vals1, trace1(startF1(b):endF1(b))];
                end
            end
            if ~isempty(vals1), avg_day1(k) = mean(vals1); end

            vals2 = [];
            for b = 1:numel(startF2)
                if endF2(b) > startF2(b)
                    vals2 = [vals2, trace2(startF2(b):endF2(b))];
                end
            end
            if ~isempty(vals2), avg_day2(k) = mean(vals2); end
        end

        % Add to pooled data
        valid = ~isnan(avg_day1) & ~isnan(avg_day2);
        pooled{c,1} = [pooled{c,1}; avg_day1(valid)];
        pooled{c,2} = [pooled{c,2}; avg_day2(valid)];
    end
end

%% Plot pooled correlations
figure('Name','Pooled Correlations','Color','w');
r_values = nan(numel(comparisons),1); % store correlations

for c = 1:numel(comparisons)
    d1 = comparisons{c}(1);
    d2 = comparisons{c}(2);

    xvals = pooled{c,1};
    yvals = pooled{c,2};

    subplot(3,7,c);
    scatter(xvals, yvals, 25, 'filled','MarkerFaceAlpha',0.5);
    xlabel(sprintf('Day %d avg activity', d1));
    ylabel(sprintf('Day %d avg activity', d2));
    title(sprintf('Day %d vs %d', d1, d2));
    axis square; grid off;

    if ~isempty(xvals)
        [r,p] = corr(xvals, yvals);
        r_values(c) = r;  % save r
        h = lsline; % regression line
        set(h,'LineWidth', 2);
        text(min(xvals), max(yvals), sprintf('r = %.2f\np = %.3g', r, p), ...
            'FontSize',12,'VerticalAlignment','top');
    end
end

%% Build correlation matrix across days
allDays = unique(days);        % list of days you’re analyzing
nDays = max(allDays);            % matrix will be nDays x nDays
corrMat = nan(nDays, nDays);     % preallocate

for c = 1:numel(comparisons)
    d1 = comparisons{c}(1);
    d2 = comparisons{c}(2);
    xvals = pooled{c,1};
    yvals = pooled{c,2};

    if ~isempty(xvals)
        r = corr(xvals, yvals);
        corrMat(d1,d2) = r;
        corrMat(d2,d1) = r; % symmetric
    end
end

for d = allDays
    corrMat(d,d) = 1; % self-correlation
end



%% Plot correlation matrix
figure('Name','Day Correlation Matrix','Color','w');
imagesc(corrMat, [0 1]);   % keep symmetric color scale
colormap(flipud(copper));
cb = colorbar;
cb.Ticks = [0, 1];
cb.TickLabels = {num2str(0), num2str(1)};
clim([0 1])
axis square;
set(gca,'XTick',allDays,'YTick',allDays,'FontSize',30, 'TickLength',[0 0]);
ylabel(cb, 'correlation (r)', 'FontSize', 30);
xlabel('day'); ylabel('day');
%exportgraphics(gca,'SevenDay_CorrMatrix.eps','ContentType','vector')

%% --- Distance Dependence Plot ---
day7 = 7;
maxLag = day7 - 1; 
avgCorrByLag = nan(1,maxLag);
semCorrByLag = nan(1,maxLag);

for lag = 1:maxLag
    vals = [];
    for d1 = 1:(day7-lag)
        d2 = d1 + lag;
        if ~isnan(corrMat(d1,d2))
            vals(end+1) = corrMat(d1,d2); %#ok<AGROW>
        end
    end
    if ~isempty(vals)
        avgCorrByLag(lag) = mean(vals);
        semCorrByLag(lag) = std(vals) / sqrt(numel(vals)); % SEM
    end
end


figure('Name','Distance Dependence','Color','w');

% plot error bars only (no marker)
h = errorbar(1:maxLag, avgCorrByLag, semCorrByLag, ...
    'LineStyle','-', 'Color',[0.55 0.35 0.05], 'LineWidth',3, ...
    'Marker','none');

hold on;

% overlay clean circle markers
plot(1:maxLag, avgCorrByLag, 'o', ...
    'MarkerSize',8, ...
    'MarkerFaceColor',[0.55 0.35 0.05], ...
    'MarkerEdgeColor',[0.55 0.35 0.05]);
xlabel('distance between days');
ylabel('mean correlation (r)');
xlim([0.5 maxLag+0.5]); ylim([0 0.7]);
set(gca,'LineWidth', 2,'XTick',allDays,'YTick',[0 0.7],'tickdir', 'out','FontSize',30);
grid off; axis square; box off;
%exportgraphics(gca,'SevenDay_Day7Distance.eps','ContentType','vector')

%% --- Spatial Footprints ---


load("/Users/nitishpatel/Desktop/Miniscope Feeding/Seven Days/D_1-2.mat"); % Change to D_1-7.mat for that comparison
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Seven Days/Mouse_D/Day 1/footprints.mat");
A1 = permute(A, [2, 3, 1]);  %  Mouse D from Day 1 
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Seven Days/Mouse_D/Day 2/footprints.mat");
A2 = permute(A, [2, 3, 1]);  %  Mouse D from Day 2
%load("/Users/nitishpatel/Desktop/Miniscope Feeding/Seven Days/Mouse_D/Day 7/footprints.mat");
%A2 = permute(A, [2, 3, 1]);  %  Mouse D from Day 7


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

% === Parameters ===
% Number of colors you need
n_colors = numel(co_registered_idx);  % number of co-registered cells or colors needed

hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb


cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);

% === Helper function ===
normalize = @(x) x / max(x(:));


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
exportgraphics(gca,'MouseD_1-2_Day1Footprints.eps','ContentType','vector')


% === Parameters ===
sigma = 2;
threshold = 0.2;
alpha_val = 0.6;
[H, W, ~] = size(A2);

% === CellReg Mapping ===
map = cell_registered_struct.cell_to_index_map;
co_registered_idx = find(map(:,1) > 0 & map(:,2) > 0);  % co-registered only
n_coreg = numel(co_registered_idx);
% === Pastel spectrum: red → violet in soft tones ===
% === Pastel spectrum: red → violet ===
hues = linspace(0, 0.83, n_coreg);   % red to violet (HSV hue)
saturation = 0.4;                    % reduced saturation = pastel
value = 1.0;                         % max brightness allowed by hsv2rgb


cmap = hsv2rgb([hues(:), repmat(saturation, n_coreg, 1), repmat(value, n_coreg, 1)]);
% Convert HSV to pastel RGB

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
exportgraphics(gca,'MouseD_1-2_Day2Footprints.eps','ContentType','vector')


