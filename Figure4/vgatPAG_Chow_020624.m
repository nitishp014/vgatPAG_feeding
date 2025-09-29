% l/vlPAG Chow Consumption
% 2/6/2024


cd('/Users/nitishpatel/Documents/MATLAB/');

ch_on_off = readmatrix('vgat_wnt_test1.xlsx','Sheet','chow','Range','M4:M13');
yfp_on_off = readmatrix('vgat_wnt_test1.xlsx','Sheet','chow','Range','M23:M32');


p = ranksum(ch_on_off, yfp_on_off); % p-value is 0.08

b = bar(0.35, [mean(ch_on_off)], 'FaceColor', 'flat');
hold on;
c = bar(0.15, [mean(yfp_on_off)], 'FaceColor', 'flat');
b.CData = [234/255 144/255 141/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [170/255 170/255 170/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.15 0.35], [mean(yfp_on_off), mean(ch_on_off)], [std(yfp_on_off)/sqrt(length(yfp_on_off)), std(ch_on_off)/sqrt(length(ch_on_off))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
plot(0.35, ch_on_off, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.15, yfp_on_off, 'ok', 'MarkerFaceColor', 'k');
hold on;
xlim([0 0.5])
yticks([-20 30])
ylim([-20 30])
legend({'l/vlPAG-Vgat ChR2', 'l/vlPAG-Vgat YFP'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'\Delta chow weight (mg)'; '(ON-OFF)'}, 'FontSize', 24);
text(0.24, 28, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 20); 
text(0.24, 24.7, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20);
text(0.24, 25.5, '_____', 'HorizontalAlignment', 'center', 'FontSize', 40);
