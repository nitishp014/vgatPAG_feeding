% l/vlPAG Devaluing Chocolate vs Cheese Day 1 Consumption (Cheese)
% 05/24/2024


cd('/Users/nitishpatel/Documents/MATLAB/');

ch_on_offS = readmatrix('vgat_wnt_test.xlsx','Sheet','cheese','Range','M123:M130');
yfp_on_offS = readmatrix('vgat_wnt_test.xlsx','Sheet','cheese','Range','M142:M149');


p = ranksum(ch_on_offS, yfp_on_offS);



b = bar(0.35, [mean(ch_on_offS)], 'FaceColor', 'flat');
hold on;
c = bar(0.15, [mean(yfp_on_offS)], 'FaceColor', 'flat');
b.CData = [234/255 144/255 141/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [170/255 170/255 170/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.15 0.35], [mean(yfp_on_offS), mean(ch_on_offS)], [std(yfp_on_offS)/sqrt(length(yfp_on_offS)), std(ch_on_offS)/sqrt(length(ch_on_offS))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 24)
plot(0.35, ch_on_offS, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.15, yfp_on_offS, 'ok', 'MarkerFaceColor', 'k');
hold on;
xlim([0 0.5])
yticks([-200 450])
ylim([-200 450])
legend({'PAG GABA ChR2-YFP', 'PAG GABA YFP'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'\Delta cheese weight (mg)'; '(ON-OFF)'}, 'FontSize', 24);
text(0.24, 465, '****', 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 499, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 480, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

