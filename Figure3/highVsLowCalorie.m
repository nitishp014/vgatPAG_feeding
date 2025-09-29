% l/vlPAG Chocolate Consumption
% 03/01/2024


cd('/Users/nitishpatel/Documents/MATLAB/');

ch_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','caloric v noncaloric','Range','E3:E11');
yfp_on_off = readmatrix('vgat_wnt_test.xlsx','Sheet','caloric v noncaloric','Range','J3:J10');

p = ranksum(ch_on_off, yfp_on_off);

figure;
b = bar(0.35, [mean(ch_on_off)], 'FaceColor', 'flat');
hold on;
c = bar(0.15, [mean(yfp_on_off)], 'FaceColor', 'flat');
b.CData = [234/255 144/255 141/255];
b.EdgeColor = 'none';
b.BarWidth = 0.15;
hold on;
c.CData = [216/255 187/255 192/255]; 
c.EdgeColor = 'none';
c.BarWidth = 0.15;
hold on;
errorbar([0.15 0.35], [mean(yfp_on_off), mean(ch_on_off)], [std(yfp_on_off)/sqrt(length(yfp_on_off)), std(ch_on_off)/sqrt(length(ch_on_off))], 'k.', 'LineWidth', 1);
hold on;
set(gca,'linewidth',2, 'TickDir','out', 'xticklabel',[], 'xtick', [], 'XAxisLocation', 'origin', 'fontsize', 30)
plot(0.35, ch_on_off, 'ok', 'MarkerFaceColor', 'k');
hold on;
plot(0.15, yfp_on_off, 'ok', 'MarkerFaceColor', 'k');
hold on;
xlim([0 0.5])
yticks([-100 250])
ylim([-100 250])
legend({'l/vlPAG vgat ChR2-YFP', '/vlPAG vgat YFP'}, 'Location', 'eastoutside');
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


ylabel({'\Delta food weight (mg)'; '(ON-OFF)'}, 'FontSize', 30);
text(0.24, 265, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
%text(0.24, 309, ['p = ' num2str(p, 2)], 'HorizontalAlignment', 'center', 'FontSize', 20); % n = 10 each groupo
text(0.24, 280, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);

