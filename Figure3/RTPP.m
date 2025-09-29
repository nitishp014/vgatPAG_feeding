chr2_on = [];
chr2_off = [];
yfp_on = [];
yfp_off = [];

cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/ChR2");
chr2_on = [chr2_on readmatrix('C85B.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C85H_rerun.xls','Range', 'W16:W16')];
%chr2_on = [chr2_on readmatrix('C86E.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C86F.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C87_0_rerun.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C87_1.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C88A_0_rerun.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C88A_1.xls','Range', 'W16:W16')];

chr2_off = [chr2_off readmatrix('C85B.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C85H_rerun.xls','Range', 'AH16:AH16')];
%chr2_off = [chr2_off readmatrix('C86E.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C86F.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C87_0_rerun.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C87_1.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C88A_0_rerun.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C88A_1.xls','Range', 'AH16:AH16')];

cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/YFP");
yfp_on = [yfp_on readmatrix('Y55_0.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y55_1.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y55_2_rerun.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y56_0.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y56_1.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y58_0_rerun.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y58_1.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y61.xls','Range', 'W16:W16')];
%yfp_on = [yfp_on readmatrix('Y85.xls','Range', 'W16:W16')];

cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/YFP");
yfp_off = [yfp_off readmatrix('Y55_0.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y55_1.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y55_2_rerun.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y56_0.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y56_1.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y58_0_rerun.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y58_1.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y61.xls','Range', 'AH16:AH16')];
%yfp_off = [yfp_off readmatrix('Y85.xls','Range', 'AH16:AH16')];



yfp_diff = yfp_on - yfp_off;
chr2_diff = chr2_on - chr2_off;
p = ranksum(yfp_diff, chr2_diff);
figure;
b = bar(0.35, [mean(chr2_diff)], 'FaceColor', 'flat');
hold on;
c = bar(0.15, [mean(yfp_diff)], 'FaceColor', 'flat');
b.CData = [234/255 144/255 141/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [170/255 170/255 170/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.15 0.35], [mean(yfp_diff), mean(chr2_diff)], [std(yfp_diff)/sqrt(length(yfp_diff)), std(chr2_diff)/sqrt(length(chr2_diff))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 30)
plot(0.35, chr2_diff, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.15, yfp_diff, 'ok', 'MarkerFaceColor', 'k');
hold on;
xlim([0 0.5])
yticks([-150 200])
ylim([-150 200])
legend({'l/vlPAG vgat ChR2-YFP', 'l/vlPAG vgat YFP'}, 'Location', 'eastoutside');
legend boxoff
box off

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

ylabel({'\Delta time stimulated side'; 'of RTPP (s) (ON-OFF)'}, 'FontSize', 30);
text(0.24, 212, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 222, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);







chr2_on = [];
chr2_off = [];
yfp_on = [];
yfp_off = [];



% Baseline Day (laser OFF on both sides, to make sure no preference on
% either side)
cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/ChR2_Baseline");
chr2_on = [chr2_on readmatrix('C85B_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C85H_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C86E_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C86F_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C87C_0_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C87C_1_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C88A_0_baseline.xls','Range', 'W16:W16')];
chr2_on = [chr2_on readmatrix('C88A_1_baseline.xls','Range', 'W16:W16')];

chr2_off = [chr2_off readmatrix('C85B_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C85H_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C86E_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C86F_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C87C_0_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C87C_1_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C88A_0_baseline.xls','Range', 'AH16:AH16')];
chr2_off = [chr2_off readmatrix('C88A_1_baseline.xls','Range', 'AH16:AH16')];

cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/YFP_Baseline");
yfp_on = [yfp_on readmatrix('Y55_0_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y55_1_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y55_2_baseline.xls','Range', 'W16:W16')];
%yfp_on = [yfp_on readmatrix('Y56_0_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y56_1_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y58_0_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y58_1_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y61_baseline.xls','Range', 'W16:W16')];
yfp_on = [yfp_on readmatrix('Y85_baseline.xls','Range', 'W16:W16')];

cd("/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/YFP_Baseline");
yfp_off = [yfp_off readmatrix('Y55_0_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y55_1_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y55_2_baseline.xls','Range', 'AH16:AH16')];
%yfp_off = [yfp_off readmatrix('Y56_0_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y56_1_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y58_0_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y58_1_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y61_baseline.xls','Range', 'AH16:AH16')];
yfp_off = [yfp_off readmatrix('Y85_baseline.xls','Range', 'AH16:AH16')];



% Example heatmap from ChR2 mouse (C85B)
% Load Excel file
filePath = '/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/C85B_position.xls';
T = readtable(filePath);

% Extract X and Y columns (columns F and G)
x = T{:,6};  % 'Body X'
y = T{:,7};  % 'Body Y'

% Remove NaNs
valid_idx = ~isnan(x) & ~isnan(y);
x = x(valid_idx);
y = y(valid_idx);

% Define chamber dimensions based on data range
x_min = floor(min(x));
x_max = ceil(max(x));
y_min = floor(min(y));
y_max = ceil(max(y));

imgWidth = x_max - x_min + 1;
imgHeight = y_max - y_min + 1;

% Shift coordinates so they fit the new image
x_shifted = round(x - x_min + 1);
y_shifted = round(y - y_min + 1);

% Initialize heatmap
heatmap = zeros(imgHeight, imgWidth);

% Accumulate occupancy
for i = 1:length(x_shifted)
    xi = x_shifted(i);
    yi = y_shifted(i);
    % Clip to image bounds
    xi = min(max(xi, 1), imgWidth);
    yi = min(max(yi, 1), imgHeight);
    heatmap(yi, xi) = heatmap(yi, xi) + 1;
end

% Smooth with 2D Gaussian filter
sigma = 6;
heatmap_blurred = imgaussfilt(heatmap, sigma);

% Normalize to range [0, 1]
heatmap_normalized = heatmap_blurred / max(heatmap_blurred(:));

% Plot heatmap
figure;
imagesc(heatmap_normalized);
axis equal tight;
n1 = 200; % length of original parula portion
n2 = 56;  % length of red extension
parula_part = parula(n1);

% Blend from yellow (end of parula) to red
start_color = parula_part(end, :); % yellow-ish
end_color = [1 0 0];               % red
ramp = [linspace(start_color(1), end_color(1), n2)', ...
        linspace(start_color(2), end_color(2), n2)', ...
        linspace(start_color(3), end_color(3), n2)'];

% Combine parula with red ramp
custom_map = [parula_part; ramp];
colormap(custom_map);
colorbar;
set(gca, 'YDir', 'normal', 'FontSize', 30);

% Remove x and y axis ticks and labels
xticks([]);
yticks([]);
xlabel('');
ylabel('');

% Add colorbar with only "min" and "max" labels
c = colorbar;
c.Ticks = [0 0.6];
c.TickLabels = {'min', 'max'};
c.Label.String = '';
c.TickLength = 0; 
c.EdgeColor = 'none';
c.Box = "off";
clim([0 0.6]);  % Optional: limit range
exportgraphics(gca,'RTPP_Heatmap.eps','ContentType','vector')


%%%Option to add a baseline day heatmap
% Load Excel file
filePath = '/Users/nitishpatel/Documents/MATLAB/VGAT_PAG_RTPP/Y55_0_baseline.xls';
T = readtable(filePath);

% Extract X and Y columns (columns F and G)
x = T{:,6};  % 'Body X'
y = T{:,7};  % 'Body Y'

% Remove NaNs
valid_idx = ~isnan(x) & ~isnan(y);
x = x(valid_idx);
y = y(valid_idx);

% Define chamber dimensions based on data range
x_min = floor(min(x));
x_max = ceil(max(x));
y_min = floor(min(y));
y_max = ceil(max(y));

imgWidth = x_max - x_min + 1;
imgHeight = y_max - y_min + 1;

% Shift coordinates so they fit the new image
x_shifted = round(x - x_min + 1);
y_shifted = round(y - y_min + 1);

% Initialize heatmap
heatmap = zeros(imgHeight, imgWidth);

% Accumulate occupancy
for i = 1:length(x_shifted)
    xi = x_shifted(i);
    yi = y_shifted(i);
    % Clip to image bounds
    xi = min(max(xi, 1), imgWidth);
    yi = min(max(yi, 1), imgHeight);
    heatmap(yi, xi) = heatmap(yi, xi) + 1;
end

% Smooth with 2D Gaussian filter
sigma = 6;
heatmap_blurred = imgaussfilt(heatmap, sigma);

% Normalize to range [0, 1]
heatmap_normalized = heatmap_blurred / max(heatmap_blurred(:));

% Plot heatmap
figure;
imagesc(heatmap_normalized);
axis equal tight;
n1 = 200; % length of original parula portion
n2 = 56;  % length of red extension
parula_part = parula(n1);

% Blend from yellow (end of parula) to red
start_color = parula_part(end, :); % yellow-ish
end_color = [1 0 0];               % red
ramp = [linspace(start_color(1), end_color(1), n2)', ...
        linspace(start_color(2), end_color(2), n2)', ...
        linspace(start_color(3), end_color(3), n2)'];

% Combine parula with red ramp
custom_map = [parula_part; ramp];
colormap(custom_map);
set(gca, 'YDir', 'normal', 'FontSize', 30);

% Remove x and y axis ticks and labels
xticks([]);
yticks([]);
xlabel('');
ylabel('');

% Add colorbar with only "min" and "max" labels
c = colorbar;
c.Ticks = [0 0.6];
c.TickLabels = {'min', 'max'};
c.Label.String = '';
c.TickLength = 0; 
c.EdgeColor = 'none';
c.Box = "off";
clim([0 0.6]);  % Optional: limit range
exportgraphics(gca,'RTPP_HeatmapBaseline.eps','ContentType','vector')
