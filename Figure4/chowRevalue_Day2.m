% l/vlPAG Chow (When Cheese is Devalued) Consumption
% 2/6/2024


cd('/Users/nitishpatel/Documents/MATLAB/');

ch_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','chow','Range','M80:M87');
yfp_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','chow','Range','M98:M105');


p = ranksum(ch_on_off, yfp_on_off); % p-value is 0.003

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
yticks([-5 45])
ylim([-5 45])
legend({'PAG GABA ChR2-YFP', 'PAG GABA YFP'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'\Delta chow weight (mg)'; '(ON-OFF)'}, 'FontSize', 24);
text(0.24, 47, '**', 'HorizontalAlignment', 'center', 'FontSize', 30); 
%text(0.24, 24.7, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20);
text(0.24, 49, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
