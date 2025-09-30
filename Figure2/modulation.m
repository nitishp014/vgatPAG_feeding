% Define root directory
rootDir = '/Users/nitishpatel/Desktop/Miniscope Feeding/';
mice = {'Mouse_A', 'Mouse_B', 'Mouse_D', 'Mouse_E', 'Mouse_F', 'Mouse_G'};
sessionName = 'Session 3';
fs = 15; % Assuming 15 Hz frame rate
preFrames = 4 * fs; % 4 seconds before eating (A)
postFrames = 4 * fs; % 4 seconds after eating starts (B)
numShuffles = 1000; % 1000 random shifts for null distribution
shiftAmount = 300; % Shift by 300 frames

% Initialize arrays to store the fractions for each mouse
positiveFractions = NaN(1, length(mice)); % Fraction of positively modulated cells for each mouse
negativeFractions = NaN(1, length(mice)); % Fraction of negatively modulated cells for each mouse

% Initialize variables to accumulate counts across all mice
totalPosModulatedCount = 0;
totalNegModulatedCount = 0;
totalCellsCount = 0; % Total number of cells across all mice

for mIdx = 1:length(mice)
    mouse = mice{mIdx};
    sessionPath = fullfile(rootDir, mouse, 'Week 1', sessionName);
    
    % Find the _combined_filt.mat file
    matFiles = dir(fullfile(sessionPath, '*_combined.mat'));
    if isempty(matFiles)
        warning('No matching .mat files found for %s', mouse);
        continue;
    end
    matFile = fullfile(sessionPath, matFiles(1).name);
    
    % Load data
    data = load(matFile);
    sumC = data.sumC; % Neural activity (N x M)
    eating_vector = data.eating_vector; % Eating events (1 x M)
    
    [numCells, numFrames] = size(sumC);
    
    % Find eating bout starts (0 → {1,2,3,4})
    eatingStarts = find(diff([0 eating_vector]) > 0);
    
    % Remove bouts too close to start or end
    validBouts = eatingStarts(eatingStarts > preFrames & eatingStarts + postFrames <= numFrames);
    
    if isempty(validBouts)
        warning('No valid eating bouts found for %s in %s', mouse, sessionName);
        continue;
    end
    
    % Compute B-A for each cell
    B_A_actual = nan(numCells, 1);
    
    for cellIdx = 1:numCells
        B_vals = nan(length(validBouts), 1);
        A_vals = nan(length(validBouts), 1);

        for i = 1:length(validBouts)
            boutStart = validBouts(i);
            B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
            A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
        end

        B_A_actual(cellIdx) = mean(B_vals - A_vals, 'omitnan');
    end
    
    % Generate null distribution using shuffled eating_vector
B_A_shuffled = nan(numCells, numShuffles);

for shufIdx = 1:numShuffles
    % Shuffle the eating_vector for each shuffle independently
    shuffled_eating = circshift(eating_vector, shiftAmount); % Random shift
    
    % Find the eating bout starts (0 → {1,2,3,4})
    shuffled_starts = find(diff([0 shuffled_eating]) > 0);
    
    % Remove bouts too close to start or end
    shuffled_validBouts = shuffled_starts(shuffled_starts > preFrames & shuffled_starts + postFrames <= numFrames);
    
    if isempty(shuffled_validBouts)
        continue;
    end
    
    for cellIdx = 1:numCells
        B_vals = nan(length(shuffled_validBouts), 1);
        A_vals = nan(length(shuffled_validBouts), 1);
        
        for i = 1:length(shuffled_validBouts)
            boutStart = shuffled_validBouts(i);
            B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
            A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
        end
        
        B_A_shuffled(cellIdx, shufIdx) = mean(B_vals - A_vals, 'omitnan');
    end
    eating_vector = shuffled_eating;
end

    
    % Compute significance thresholds
    perc10 = prctile(B_A_shuffled, 15, 2);
    perc90 = prctile(B_A_shuffled, 85, 2);
    
    % Count significantly modulated cells for each mouse
    posModulatedCount = sum(B_A_actual >= perc90, 'omitnan');
    negModulatedCount = sum(B_A_actual <= perc10, 'omitnan');
   
    % Accumulate counts
    totalPosModulatedCount = totalPosModulatedCount + posModulatedCount;
    totalNegModulatedCount = totalNegModulatedCount + negModulatedCount;
    totalCellsCount = totalCellsCount + numCells;

    % Calculate the fraction for each mouse
    positiveFractions(mIdx) = posModulatedCount / numCells;
    negativeFractions(mIdx) = negModulatedCount / numCells;
end

% Calculate the average fractions across all mice
avgPositiveFraction = mean(positiveFractions, 'omitnan');
avgNegativeFraction = mean(negativeFractions, 'omitnan');


% Compute non-modulated fraction
avgNonModulatedFraction = 1 - (avgPositiveFraction + avgNegativeFraction);

% Create labels and data
labels = {'Negative - 37%', 'Positive - 23%', 'Non-modulated - 40%'};
fractions = [avgNegativeFraction, avgPositiveFraction, avgNonModulatedFraction];

% Create pie chart
figure;
pie(fractions, labels);
%exportgraphics(gca,'ModulatedCells_PieChart.eps','ContentType','vector')



% Calculate SEM for each group
sem_neg = std(negativeFractions) / sqrt(length(negativeFractions));
sem_pos = std(positiveFractions) / sqrt(length(positiveFractions));


% Plot the average fraction of positively and negatively modulated cells
figure('Position',[100, 100, 375, 425]);
hold on;
bar(0.15, avgNegativeFraction, 'FaceColor', [190, 211, 237] / 255,'EdgeColor','none', 'BarWidth', 0.15); % Bar for negative modulation
bar(0.35, avgPositiveFraction, 'FaceColor', [106, 145, 247] / 255,'EdgeColor','none', 'BarWidth', 0.15); % Bar for positive modulation

% Plot SEM as error bars
errorbar(0.15, mean(negativeFractions), sem_neg, 'k.', 'LineWidth', 1);
errorbar(0.35, mean(positiveFractions), sem_pos, 'k.', 'LineWidth', 1);

% Plot individual data points
plot(0.15, negativeFractions, 'ok', 'MarkerFaceColor', 'k');
plot(0.35, positiveFractions, 'ok', 'MarkerFaceColor', 'k');




% Formatting the plot
xlim([0 0.5]);
ylim([0, 0.5])
yticks([0, 0.5])
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
ylabel('fraction all cells');
legend({'eat (-)', 'eat (+)'}, 'Location', 'eastoutside');
legend boxoff
box off
hold off;

% Perform statistical test
p = ranksum(negativeFractions, positiveFractions);

% Add asterisks for significance
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

text(0.24, 0.52, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
text(0.24, 0.54, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
%exportgraphics(gca,'negvspos_quantification.eps','ContentType','vector')

% Display the total counts
disp(['Total positively modulated cells: ', num2str(totalPosModulatedCount)]);
disp(['Total negatively modulated cells: ', num2str(totalNegModulatedCount)]);
disp(['Total cells across all mice: ', num2str(totalCellsCount)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Below is best examples of negatively modulated mice from each session
   %G12 - 45, 
   %G13 - 21
   %F12 - 31, 37, 49!, 53, 
   %F13 - 11, 54, 
   %E12 - 14, 
   %E13 - 14, 
   %D12 - 39, 54, 77
   %B12 - 40
   %B13 - 9, 
   %A13 - 13 

 % Below is best examples of positively modulated mice from each session
   %G13 - 2, 6, 26, 28, 32
   %F13 - 72,73
   %E13 - 5, 8, 25, 
   %D13 - 17!, 55, 82
   %B13 - 28

%nonmod cells?
%G13 - 12, 
   
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_F/Week 1/Session 2/F12_combined_filt.mat")
sumC = sumC;
eating_vector = eating_vector;
cell_idx = 49;
sampling_rate = 15;
trace_color = [95 83 130] / 255;
sem_shade_color = [224 216 238] / 255;

  % Ensure eating_vector is a row vector
    eating_vector = eating_vector(:)'; 

    % Find eating bout start times (when 0 -> 1)
    eating_starts = find(diff(eating_vector > 0) == 1); 

    % Define time window (5 seconds before and after eating start)
    window_size = 5 * sampling_rate; % Convert seconds to frames
    
    % Initialize matrix to store valid bout-aligned activity
    valid_bouts = [];

    % Extract activity around eating starts
    for i = 1:length(eating_starts)
        start_idx = eating_starts(i) - window_size;
        end_idx = eating_starts(i) + window_size;

        % Ensure indices are within valid bounds
        if start_idx > 0 && end_idx <= size(sumC, 2)
            valid_bouts = [valid_bouts; sumC(cell_idx, start_idx:end_idx)];
        end
    end

    % Compute mean and SEM
    avg_activity = mean(valid_bouts, 1);
    sem_activity = std(valid_bouts, [], 1) ./ sqrt(size(valid_bouts, 1));




    % Time axis (5 seconds before to 5 seconds after)
    time_axis = linspace(-5, 5, length(avg_activity));

    % Plot the mean activity
    figure; hold on;
        % Create shaded error region
    fill([time_axis, fliplr(time_axis)], ...
         [avg_activity - sem_activity, fliplr(avg_activity + sem_activity)], ...
         sem_shade_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName','SEM');
    plot(time_axis, avg_activity, 'k', 'LineWidth', 2, 'Color', trace_color, 'DisplayName','Avg Trace');



        
        % Add a vertical red solid line at 0 seconds (start of behavior)
        plot([0, 0], [-0.1, 3], 'r--', 'LineWidth', 2, 'DisplayName', 'Start of Behavior');
        
        % Make the horizontal 0 line lighter (gray instead of black)
        plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, 'DisplayName', '0 df/F');

        
        % Darken the x and y axis lines
        ax = gca;
        ax.XAxis.LineWidth = 2;  % Darker x-axis line
        ax.YAxis.LineWidth = 2;  % Darker y-axis line

        % Set axis limits and ticks
        xlim([-5, 5]);  % 10 seconds before and after behavior onset
        ylim([-0.1, 3]);  % Y-axis range
        xticks([-5, 0, 5]);  % Only label -5, 0, and 5 on the x-axis
        yticks([-0.1, 3]);  % Y-axis labels
        yticklabels({'-0.1', '3'});
        
        % Set additional plot aesthetics
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);
        
        % Add labels
        xlabel('eating onset (s)');
        ylabel('mean df/F (z-scored)');
        title(sprintf('Example Cell #%d', cell_idx), 'FontSize', 15);
        % Finalize plot
        legend off;
        hold off;
        exportgraphics(gca,'negcell_example.eps','ContentType','vector')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   % Below is best examples of negatively modulated mice from each session
   %G12 - 45, 
   %G13 - 21
   %F12 - 31, 37, 49!, 53, 
   %F13 - 11, 54, 
   %E12 - 14, 
   %E13 - 14, 
   %D12 - 39, 54, 77
   %B12 - 40
   %B13 - 9, 
   %A13 - 13 

 % Below is best examples of positively modulated mice from each session
   %G13 - 2, 6, 26, 28, 32
   %F13 - 72,73
   %E13 - 5, 8, 25, 
   %D13 - 17!, 55, 82
   %B13 - 28

%nonmod cells?
%G13 - 12, 
   
load("/Users/nitishpatel/Desktop/Miniscope Feeding/Mouse_D/Week 1/Session 3/D13_combined_filt.mat")
sumC = sumC;
eating_vector = eating_vector;
cell_idx = 17;
sampling_rate = 15;
trace_color = [95 83 130] / 255;
sem_shade_color = [224 216 238] / 255;

  % Ensure eating_vector is a row vector
    eating_vector = eating_vector(:)'; 

    % Find eating bout start times (when 0 -> 1)
    eating_starts = find(diff(eating_vector > 0) == 1); 

    % Define time window (5 seconds before and after eating start)
    window_size = 5 * sampling_rate; % Convert seconds to frames
    
    % Initialize matrix to store valid bout-aligned activity
    valid_bouts = [];

    % Extract activity around eating starts
    for i = 1:length(eating_starts)
        start_idx = eating_starts(i) - window_size;
        end_idx = eating_starts(i) + window_size;

        % Ensure indices are within valid bounds
        if start_idx > 0 && end_idx <= size(sumC, 2)
            valid_bouts = [valid_bouts; sumC(cell_idx, start_idx:end_idx)];
        end
    end

    % Compute mean and SEM
    avg_activity = mean(valid_bouts, 1);
    sem_activity = std(valid_bouts, [], 1) ./ sqrt(size(valid_bouts, 1));




    % Time axis (5 seconds before to 5 seconds after)
    time_axis = linspace(-5, 5, length(avg_activity));

    % Plot the mean activity
    figure; hold on;
        % Create shaded error region
    fill([time_axis, fliplr(time_axis)], ...
         [avg_activity - sem_activity, fliplr(avg_activity + sem_activity)], ...
         sem_shade_color, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName','SEM');
    plot(time_axis, avg_activity, 'k', 'LineWidth', 2, 'Color', trace_color, 'DisplayName','Avg Trace');



        
        % Add a vertical red solid line at 0 seconds (start of behavior)
        plot([0, 0], [-0.1, 3], 'r--', 'LineWidth', 2, 'DisplayName', 'Start of Behavior');
        
        % Make the horizontal 0 line lighter (gray instead of black)
        plot(xlim, [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, 'DisplayName', '0 df/F');

        
        % Darken the x and y axis lines
        ax = gca;
        ax.XAxis.LineWidth = 2;  % Darker x-axis line
        ax.YAxis.LineWidth = 2;  % Darker y-axis line

        % Set axis limits and ticks
        %xlim([-5, 5]);  % 10 seconds before and after behavior onset
        %ylim([-2, 3]);  % Y-axis range
        %xticks([-5, 0, 5]);  % Only label -5, 0, and 5 on the x-axis
        %yticks([-2, 3]);  % Y-axis labels
        %yticklabels({'-2', '3'});
        
        % Set additional plot aesthetics
        set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 24);
        
        % Add labels
        xlabel('eating onset (s)');
        ylabel('mean df/F (z-scored)');
        title(sprintf('Example Cell #%d', cell_idx), 'FontSize', 15);
        % Finalize plot
        legend off;
        hold off;
       % exportgraphics(gca,'negcell_example.eps','ContentType','vector')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Proportion of Positively and Negatively Modulated Cells by Food Type -
% Likely not very useful because of small sample sizes (few numbers of
% eating bouts per food)

% Define root directory
rootDir = '/Users/nitishpatel/Desktop/Miniscope Feeding/';
mice = {'Mouse_A', 'Mouse_B', 'Mouse_D', 'Mouse_E', 'Mouse_F', 'Mouse_G'};
sessionName = 'Session 3';
fs = 15; % Assuming 15 Hz frame rate
preFrames = 4 * fs; % 4 seconds before eating (A)
postFrames = 4 * fs; % 4 seconds after eating starts (B)
numShuffles = 1000; % 300 random shifts for null distribution
shiftAmount = 300; % Shift by 300 frames

% Food types
foodTypes = {'Cheese', 'Broccoli', 'Chocolate', 'Carrot'};
numFoodTypes = length(foodTypes);

% Initialize arrays to store the fractions for each mouse and food type
positiveFractions = NaN(length(mice), numFoodTypes);
negativeFractions = NaN(length(mice), numFoodTypes);

% Initialize variables to accumulate counts across all mice and food types
totalPosModulatedCount = zeros(1, numFoodTypes);
totalNegModulatedCount = zeros(1, numFoodTypes);
totalCellsCount = zeros(1, numFoodTypes);

for mIdx = 1:length(mice)
    mouse = mice{mIdx};
    sessionPath = fullfile(rootDir, mouse, 'Week 1', sessionName);
    
    % Find the _combined_filt.mat file
    matFiles = dir(fullfile(sessionPath, '*_combined_filt.mat'));
    if isempty(matFiles)
        warning('No matching .mat files found for %s', mouse);
        continue;
    end
    matFile = fullfile(sessionPath, matFiles(1).name);
    
    % Load data
    data = load(matFile);
    sumC = data.sumC; % Neural activity (N x M)
    eating_vector = data.eating_vector; % Eating events (1 x M)
    
    [numCells, numFrames] = size(sumC);
    
    % Iterate over each food type
    for fIdx = 1:numFoodTypes
        foodCode = fIdx; % Food codes are 1-based
        foodName = foodTypes{fIdx};
        
        % Find eating bout starts for the specific food type
        eatingStarts = find(diff([0 eating_vector == foodCode]) == 1);
        
        % Remove bouts too close to start or end
        validBouts = eatingStarts(eatingStarts > preFrames & eatingStarts + postFrames <= numFrames);
        
        if isempty(validBouts)
            warning('No valid eating bouts found for %s in %s for food type %s', mouse, sessionName, foodName);
            continue;
        end
        
        % Compute B-A for each cell
        B_A_actual = nan(numCells, 1);
        
        for cellIdx = 1:numCells
            B_vals = nan(length(validBouts), 1);
            A_vals = nan(length(validBouts), 1);
            
            for i = 1:length(validBouts)
                boutStart = validBouts(i);
                B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
                A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
            end
            
            B_A_actual(cellIdx) = mean(B_vals - A_vals, 'omitnan');
        end
        
        % Generate null distribution using shuffled eating_vector
        B_A_shuffled = nan(numCells, numShuffles);
        
        for shufIdx = 1:numShuffles
            % Shuffle the eating_vector for each shuffle independently
            shuffled_eating = circshift(eating_vector, shiftAmount); 
            
            % Find the eating bout starts for the specific food type
            shuffled_starts = find(diff([0 shuffled_eating == foodCode]) == 1);
            
            % Remove bouts too close to start or end
            shuffled_validBouts = shuffled_starts(shuffled_starts > preFrames & shuffled_starts + postFrames <= numFrames);
            
            if isempty(shuffled_validBouts)
                continue;
            end
            
            for cellIdx = 1:numCells
                B_vals = nan(length(shuffled_validBouts), 1);
                A_vals = nan(length(shuffled_validBouts), 1);
                
                for i = 1:length(shuffled_validBouts)
                    boutStart = shuffled_validBouts(i);
                    B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
                    A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
                end
                
                B_A_shuffled(cellIdx, shufIdx) = mean(B_vals - A_vals, 'omitnan');
            end
            eating_vector = shuffled_eating;
        end
        
        % Compute significance thresholds
        perc10 = prctile(B_A_shuffled, 20, 2);
        perc90 = prctile(B_A_shuffled, 80, 2);
        
        % Count significantly modulated cells for each mouse and food type
        posModulatedCount = sum(B_A_actual >= perc90, 'omitnan');
        negModulatedCount = sum(B_A_actual <= perc10, 'omitnan');
        
        % Accumulate counts
        totalPosModulatedCount(fIdx) = totalPosModulatedCount(fIdx) + posModulatedCount;
        totalNegModulatedCount(fIdx) = totalNegModulatedCount(fIdx) + negModulatedCount;
        totalCellsCount(fIdx) = totalCellsCount(fIdx) + numCells;
        
        % Calculate the fraction for each mouse and food type
        positiveFractions(mIdx, fIdx) = posModulatedCount / numCells;
        negativeFractions(mIdx, fIdx) = negModulatedCount / numCells;
    end
end

% Calculate the average fractions across all mice for each food type
avgPositiveFraction = mean(positiveFractions, 1, 'omitnan');
avgNegativeFraction = mean(negativeFractions, 1, 'omitnan');

% Calculate SEM for each group
sem_neg = std(negativeFractions, 0, 1, 'omitnan') / sqrt(length(mice));
sem_pos = std(positiveFractions, 0, 1, 'omitnan') / sqrt(length(mice));

% Plot the average fraction of positively and negatively modulated cells by food type
figure('Position', [100, 100, 800, 600]);
hold on;
barWidth = 0.35;
groupWidth = min(0.8, barWidth * numFoodTypes);
ylim([0 0.5]);
yticks([0 0.5]);
for fIdx = 1:numFoodTypes
    x = fIdx;
    bar(x - 0.15, avgNegativeFraction(fIdx), 'FaceColor', [190, 211, 237] / 255, 'BarWidth', 0.3);
    bar(x + 0.15, avgPositiveFraction(fIdx), 'FaceColor', [106, 145, 247] / 255, 'BarWidth', 0.3);
    
    % Error bars
    errorbar(x - 0.15, avgNegativeFraction(fIdx), sem_neg(fIdx), 'k.', 'LineWidth', 1);
    errorbar(x + 0.15, avgPositiveFraction(fIdx), sem_pos(fIdx), 'k.', 'LineWidth', 1);
end

% Formatting
xticks(1:numFoodTypes);
xticklabels(foodTypes);
ylabel('Fraction of all cells');
legend({'Negative Modulation', 'Positive Modulation'}, 'Location', 'northeastoutside');
legend boxoff;
box off;
set(gca, 'FontSize', 14, 'LineWidth', 2);
hold off;

% Display the total counts
for fIdx = 1:numFoodTypes
    fprintf('%s: Positively Modulated = %d, Negatively Modulated = %d, Total Cells = %d\n', ...
        foodTypes{fIdx}, totalPosModulatedCount(fIdx), totalNegModulatedCount(fIdx), totalCellsCount(fIdx));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculates average number of eating bouts per type of food
% Define base directory
baseDir = "/Users/nitishpatel/Desktop/Miniscope Feeding/";

% List of mouse folders
mouseFolders = ["Mouse_A", "Mouse_B", "Mouse_D", "Mouse_E", "Mouse_F", "Mouse_G"];

% Session to analyze
sessionName = "Session 1"; % Change this for different sessions

% Initialize storage for bout counts
numMice = numel(mouseFolders);
boutCounts = zeros(numMice, 4); % Columns: [cheese, broccoli, chocolate, carrot]

% Loop through each mouse
for m = 1:numMice
    % Construct file path
    dataFile = fullfile(baseDir, mouseFolders(m), "Week 1", sessionName, "*_combined_filt.mat");
    matFiles = dir(dataFile);
    
    % Check if file exists
    if isempty(matFiles)
        warning("No matching file found for %s", mouseFolders(m));
        continue;
    end
    
    % Load data
    data = load(fullfile(matFiles.folder, matFiles.name), 'eating_vector');
    eating_vector = data.eating_vector;
    
    % Process eating bouts
    for foodType = 1:4
        % Find transitions from 0 to the current food type
        transitions = find(eating_vector(1:end-1) == 0 & eating_vector(2:end) == foodType);
        boutCounts(m, foodType) = numel(transitions);
    end
end

% Compute average number of bouts per food type
avgBoutCounts = mean(boutCounts, 1, 'omitnan');

% Compute SEM (Standard Error of the Mean)
semBoutCounts = std(boutCounts, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(boutCounts), 1));

% Plot settings
foodLabels = ["Cheese", "Broccoli", "Chocolate", "Carrot"];
colors = [253/255 245/255 217/255;  % Cheese 
          214/255 254/255 208/255;    % Broccoli 
          231/255 215/255 196/255; % Chocolate
          252/255 236/255 207/255];   % Carrot 
xVals = 1:4;

% Create bar graph with error bars
figure('Position',[100, 100, 375, 425]);; hold on;
barHandles = bar(xVals, avgBoutCounts, 'FaceColor', 'flat', 'EdgeColor','none');

% Apply colors to bars
for i = 1:4
    barHandles.CData(i, :) = colors(i, :);
end

% Plot error bars
errorbar(xVals, avgBoutCounts, semBoutCounts, 'k.', 'LineWidth', 1);

% Overlay individual data points (each mouse) with jitter
jitterAmount = 0.15; % Amount of horizontal spread for visibility
for i = 1:4
    xJitter = xVals(i) + (rand(numMice, 1) - 0.5) * jitterAmount; % Small random shift
    scatter(xJitter, boutCounts(:,i), 'MarkerFaceColor', 'k', 'MarkerEdgeColor','none');
end

% Labels & formatting
ylim([0 10]);
yticks([0 10]);
xticks(xVals);
xticklabels(foodLabels);
ylabel('Number of Bouts');
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 18);

% Legend
box off;
legend off;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates average duration of a bout per type of food

% Define base directory
baseDir = "/Users/nitishpatel/Desktop/Miniscope Feeding/";

% List of mouse folders
mouseFolders = ["Mouse_A", "Mouse_B", "Mouse_D", "Mouse_E", "Mouse_F", "Mouse_G"];

% Session to analyze
sessionName = "Session 2";

% Define frame rate (frames per second)
fps = 15; 

% Initialize storage for bout durations
numMice = numel(mouseFolders);
boutDurations = cell(numMice, 4); % Store durations for each bout

% Loop through each mouse
for m = 1:numMice
    % Construct file path
    dataFile = fullfile(baseDir, mouseFolders(m), "Week 1", sessionName, "*_combined_filt.mat");
    matFiles = dir(dataFile);
    
    % Check if file exists
    if isempty(matFiles)
        warning("No matching file found for %s", mouseFolders(m));
        continue;
    end
    
    % Load data
    data = load(fullfile(matFiles.folder, matFiles.name), 'eating_vector');
    eating_vector = data.eating_vector;
    
    % Process eating bouts
    for foodType = 1:4
        % Find transitions from 0 to the current food type (start of a bout)
        startIndices = find(eating_vector(1:end-1) == 0 & eating_vector(2:end) == foodType) + 1;
        
        % Find corresponding bout end indices
        durations = [];
        for startIdx = startIndices
            % Find where eating stops (next 0 or food switch)
            endIdx = find(eating_vector(startIdx:end) ~= foodType, 1, 'first');
            if isempty(endIdx)
                endIdx = length(eating_vector) - startIdx + 1; % If no end found, bout runs to end
            else
                endIdx = startIdx + endIdx - 2; % Adjust to actual index
            end
            durations = [durations; (endIdx - startIdx + 1) / fps]; % Convert to seconds
        end
        
        % Store bout durations
        boutDurations{m, foodType} = durations;
    end
end

% Compute average bout duration per food type per mouse (in seconds)
avgDurationsPerMouse = nan(numMice, 4);
for m = 1:numMice
    for foodType = 1:4
        if ~isempty(boutDurations{m, foodType})
            avgDurationsPerMouse(m, foodType) = mean(boutDurations{m, foodType});
        end
    end
end

% Compute group-level averages and SEM (in seconds)
avgBoutDurations = mean(avgDurationsPerMouse, 1, 'omitnan');
semBoutDurations = std(avgDurationsPerMouse, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(avgDurationsPerMouse), 1));

% Plot settings
foodLabels = ["Cheese", "Broccoli", "Chocolate", "Carrot"];
colors = [253/255 245/255 217/255;  % Cheese
          214/255 254/255 208/255; % Broccoli
          231/255 215/255 196/255; % Chocolate
          252/255 236/255 207/255]; % Carrot
xVals = 1:4;

% Create bar graph with error bars
figure('Position', [100, 100, 375, 425]); hold on;
barHandles = bar(xVals, avgBoutDurations, 'FaceColor', 'flat', 'EdgeColor', 'none');

% Apply colors to bars
for i = 1:4
    barHandles.CData(i, :) = colors(i, :);
end

% Plot error bars
errorbar(xVals, avgBoutDurations, semBoutDurations, 'k.', 'LineWidth', 1);

% Overlay individual data points (each mouse) with jitter
jitterAmount = 0.15; % Amount of horizontal spread for visibility
for i = 1:4
    xJitter = xVals(i) + (rand(numMice, 1) - 0.5) * jitterAmount; % Small random shift
    scatter(xJitter, avgDurationsPerMouse(:, i), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
end

% Labels & formatting
ylim([0, 40]);
yticks([0, 40]);
xticks(xVals);
xticklabels(foodLabels);
ylabel('Bout Duration (seconds)');
set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontSize', 18);

% Legend
box off;
legend off;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code displays average B-A values (bout - pre bout activity) for 

% Define root directory
rootDir = '/Users/nitishpatel/Desktop/Miniscope Feeding/';
mice = {'Mouse_A', 'Mouse_B', 'Mouse_D', 'Mouse_E', 'Mouse_F', 'Mouse_G'};
sessionName = 'Session 3';
fs = 15; % Frame rate (Hz)
preFrames = 4 * fs; % 4 seconds before eating (A)
postFrames = 4 * fs; % 4 seconds after eating starts (B)
numShuffles = 1000; % Number of random shifts for null distribution
shiftAmount = 300; % Shift amount in frames

% Food types
foodTypes = {'Cheese', 'Broccoli', 'Chocolate', 'Carrot'};
numFoodTypes = length(foodTypes);

% Initialize storage for negatively modulated cells
negModulatedBA = [];

for mIdx = 1:length(mice)
    mouse = mice{mIdx};
    sessionPath = fullfile(rootDir, mouse, 'Week 1', sessionName);
    
    % Find the _combined_filt.mat file
    matFiles = dir(fullfile(sessionPath, '*_combined_filt.mat'));
    if isempty(matFiles)
        warning('No matching .mat files found for %s', mouse);
        continue;
    end
    matFile = fullfile(sessionPath, matFiles(1).name);
    
    % Load data
    data = load(matFile);
    sumC = data.sumC; % Neural activity (N x M)
    eating_vector = data.eating_vector; % Eating events (1 x M)
    
    [numCells, numFrames] = size(sumC);
    
    % Compute B-A values for each food type
    B_A_actual = nan(numCells, numFoodTypes);
    
    for foodType = 1:numFoodTypes
        % Find eating bouts for the current food type
        eatingStarts = find(diff([0 eating_vector == foodType]) == 1);
        
        % Remove bouts too close to start or end
        validBouts = eatingStarts(eatingStarts > preFrames & eatingStarts + postFrames <= numFrames);
        
        if isempty(validBouts)
            warning('No valid eating bouts found for %s in %s for %s', mouse, sessionName, foodTypes{foodType});
            continue;
        end
        
        % Compute B-A values for each cell
        for cellIdx = 1:numCells
            B_vals = nan(length(validBouts), 1);
            A_vals = nan(length(validBouts), 1);
            
            for i = 1:length(validBouts)
                boutStart = validBouts(i);
                B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
                A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
            end
            
            B_A_actual(cellIdx, foodType) = mean(B_vals - A_vals, 'omitnan');
        end
    end
    
    % Generate null distribution using shuffled eating_vector
    B_A_shuffled = nan(numCells, numFoodTypes, numShuffles);
    
    for shufIdx = 1:numShuffles
        % Shuffle eating_vector
        shuffled_eating = circshift(eating_vector, shiftAmount);
        for foodType = 1:numFoodTypes
            shuffledStarts = find(diff([0 shuffled_eating == foodType]) == 1);
            shuffled_validBouts = shuffledStarts(shuffledStarts > preFrames & shuffledStarts + postFrames <= numFrames);
            
            if isempty(shuffled_validBouts)
                continue;
            end
            
            for cellIdx = 1:numCells
                B_vals = nan(length(shuffled_validBouts), 1);
                A_vals = nan(length(shuffled_validBouts), 1);
                
                for i = 1:length(shuffled_validBouts)
                    boutStart = shuffled_validBouts(i);
                    B_vals(i) = mean(sumC(cellIdx, boutStart:boutStart+postFrames-1), 'omitnan');
                    A_vals(i) = mean(sumC(cellIdx, boutStart-preFrames:boutStart-1), 'omitnan');
                end
                
                B_A_shuffled(cellIdx, foodType, shufIdx) = mean(B_vals - A_vals, 'omitnan');
            end
        end
        eating_vector = shuffled_eating;
    end
    
    % Compute significance thresholds
    perc10 = prctile(B_A_shuffled, 10, 3);
    %perc90 = prctile(B_A_shuffled, 90, 3); %Uncomment this for posMod
    
    % Identify negatively modulated cells
    negModulatedCells = B_A_actual <= perc10;
    %negModulatedCells = B_A_actual >= perc10;  %Uncomment this for posMod
    
    % Store B-A values for negatively modulated cells
    for cellIdx = 1:numCells
        if any(negModulatedCells(cellIdx, :))
            negModulatedBA = [negModulatedBA; B_A_actual(cellIdx, :)];
        end
    end
end

% Display the final matrix
disp('B-A matrix for negatively modulated cells (N x 4):');
disp(negModulatedBA);


