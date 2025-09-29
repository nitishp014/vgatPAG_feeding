
cd("/Users/nitishpatel/Documents/MATLAB")

y = readmatrix('vgat_wnt_test.xlsx','Sheet','food_rank','Range','S15:S70');
x = readmatrix('vgat_wnt_test.xlsx','Sheet','food_rank','Range','Q15:Q70');

food = readmatrix('vgat_wnt_test.xlsx','Sheet','food_rank','Range','X48:AE54');
categories = 1:0.5:4.5; % Categorical data on the x-axis
data = food; % Random data for 7 mice (rows) across 8 categories (columns)

% Define a colormap (e.g., using 'lines' which provides distinct colors)
%colors = lines(7); % 7 colors for 7 mice

% Create the scatter plot
figure('Position', [100, 100, 1100, 600]);  % wider figure window
hold on;
for mouse = 1:7
    scatter(categories, data(mouse, :), 90, [0.780, 0.480, 0.471], 'filled');
end
hold on;
p = polyfit(x, y, 1); % 1 indicates a first-degree polynomial (linear)

% Generate y values based on the fitted model
y_fit = polyval(p, x);

% Plot the trendline
plot(x, y_fit, 'Color',[0.780, 0.480, 0.471], 'LineWidth', 1.5); % Red line with a width of 2

% Add labels and title
ylabel({'\Delta consumption (mg)'; '(ON-OFF)'}, 'FontSize', 24);
%text(4, 400, sprintf("R^2 = %d", R_squared));
%title('Relative Ranking of Foods');
%legend('Data', 'Trendline');
%hold off;
set(gca,'linewidth',2, 'TickDir','out', 'XAxisLocation', 'origin', 'fontsize', 24)
ylim([-30 450])
yticks([0 450]);
xlim([0.7 4.8])
set(gca, 'XTick', 1:0.5:4.5); % Ensures x-ticks are at 1,2,...,8
hold on;
xticklabels({'Cherry','Carrots','Pepper','Sugar','Broccoli','Sausage','Chocolate','Cheese'});

%xticklabels({'Cherry','Carrots','Pepper','Sugar','Broccoli','Sausage','Chocolate','Cheese'});
xtickangle(70);
[r, p] = corr(x, y, 'Type', 'Spearman');

text(1.75, 350, sprintf('r = %.2f', r), 'HorizontalAlignment', 'center', 'FontSize', 20);
text(1.75, 325, sprintf('p = %.3f', p), 'HorizontalAlignment', 'center', 'FontSize', 20);


box off;
hold off;

