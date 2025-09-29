% l/vlPAG Devaluing Cheese vs Chocolate Day 2 Consumption (Cheese)
% 05/17/2024


cd('/Users/nitishpatel/Documents/MATLAB/');

ch_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','cheese','Range','M84:M91');
yfp_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','cheese','Range','M103:M110');

p = ranksum(ch_on_off, yfp_on_off);


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
yticks([-100 450])
ylim([-100 450])
legend({'PAG GABA ChR2-YFP', 'PAG GABA YFP'}, 'Location', 'eastoutside');
legend boxoff
box off

ylabel({'\Delta cheese weight (mg)'; '(ON-OFF)'}, 'FontSize', 24);
%text(0.24, 330, 'n.s.', 'HorizontalAlignment', 'center', 'FontSize', 20); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
%text(0.24, 309, '_____', 'HorizontalAlignment', 'center', 'FontSize', 40);

