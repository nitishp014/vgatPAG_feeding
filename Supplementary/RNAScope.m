% Data
ghsr = [6.5 9.3 9.4];   % % vgat cells with ghrelin receptor
lepr = [6.5 11.7 7.7];  % % vgat cells with leptin receptor

% Calculate means and SEM
means = [mean(ghsr) mean(lepr)];
sems  = [std(ghsr)/sqrt(length(ghsr)) std(lepr)/sqrt(length(lepr))];

% X positions for bars
x = [0.15 0.35];


% ---- Main panel ----
figure; hold on;

% Bars with custom colors
bar(0.35, means(1), 'FaceColor', [194 59 59]/255, 'EdgeColor','none', 'BarWidth',0.15);
hold on;
bar(0.15, means(2), 'FaceColor', [212 175 55]/255, 'EdgeColor','none', 'BarWidth',0.15);
hold on;
% Add error bars (SEM)
errorbar(x, means, sems, 'k.', 'LineWidth', 1);

plot(0.15, lepr, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.35, ghsr, 'ok', 'MarkerFaceColor', 'k');
hold on;

% Aesthetics
xlim([0 0.5]);
ylim([0 100]);  % Keep scale small to emphasize rarity
yticks([0 100]);
xticks([0.15 0.35]);
xticklabels({'{\itLepr}', '{\itGhsr}'});
ylabel('% vgat cells expressing receptor');
set(gca,'linewidth',2, 'TickDir','out', 'XAxisLocation', 'origin', 'fontsize', 24);
hold off;
box off;
%exportgraphics(gca,'rnascope_ghsrlepr.eps','ContentType','vector')