%% Load behavior
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_A/Week 2/Session 2/hungryBehavior.mat");
%% Setup
foods = {'carrot','cheese','broccoli','chocolate'};
% pastel RGB colors (normalized)
colors = [255 179 71; 255 239 130; 144 238 144; 205 133 63] ./ 255;

mice = fieldnames(behavior);
numMice = numel(mice);
numFoods = numel(foods);

% Preallocate matrices (rows = mice, cols = foods)
numBouts    = nan(numMice,numFoods);
boutLength  = nan(numMice,numFoods);
latencyFirst= nan(numMice,numFoods);

% Fill matrices safely (if a field is missing we leave NaN)
for m = 1:numMice
    for f = 1:numFoods
        % numBouts
        if isfield(behavior.(mice{m}),'numBouts') && isfield(behavior.(mice{m}).numBouts, foods{f})
            numBouts(m,f) = behavior.(mice{m}).numBouts.(foods{f});
        end
        % boutLength
        if isfield(behavior.(mice{m}),'boutLength') && isfield(behavior.(mice{m}).boutLength, foods{f})
            boutLength(m,f) = behavior.(mice{m}).boutLength.(foods{f});
        end
        % latencyFirst
        if isfield(behavior.(mice{m}),'latencyFirst') && isfield(behavior.(mice{m}).latencyFirst, foods{f})
            latencyFirst(m,f) = behavior.(mice{m}).latencyFirst.(foods{f});
        end
    end
end

%% Plot: 1x3 subplots (Number of bouts | Bout duration | Latency)
%% Measures and labels
measures = {numBouts, boutLength, latencyFirst};
ylabels  = {'number of bouts','bout duration (s)','latency to first bout (s)'};

% y-axis limits + ticks for each measure
ylims = {[0 8], [0 180], [0 200]};
yticks_all = {[0 8], [0 180], [0 200]};

for k = 1:3
    data = measures{k};               % mice x foods
    numFoods = size(data,2);          % detect number of foods
    
    % Compute mean and SEM across mice
    means = mean(data,1,'omitnan');   % 1 x foods
    sems  = std(data,0,1,'omitnan') ./ sqrt(sum(~isnan(data),1));
    
    % Create figure
    figure('Color','w','Position',[400 300 600 500],'Visible','on'); 
    hold on; box off;
    
    % Bar plot
    hb = bar(1:numFoods, means, 0.6, 'FaceColor', 'flat','EdgeColor','none');
    hb.CData = colors(1:numFoods,:);   % ensure colors match number of foods
    
    % SEM errorbars
    errorbar(1:numFoods, means, sems, 'k.','LineWidth',1);
    
    % Scatter individual mouse points (no jitter, black)
    for f = 1:numFoods
        y = data(:,f);
        idx = ~isnan(y);
        if any(idx)
            x = repmat(f,sum(idx),1);  % exact x = bar center
            scatter(x, y(idx), 60, 'k','filled');  % black dots
        end
    end
    
    % Axis formatting
    xlim([0.5 numFoods+0.5]);
    ylim(ylims{k});
    set(gca,'YTick',yticks_all{k});
    set(gca,'XTick',1:numFoods,'XTickLabel',foods,'FontSize',24,'LineWidth',2);
    ylabel(ylabels{k}, 'FontSize',24);
end

drawnow;
