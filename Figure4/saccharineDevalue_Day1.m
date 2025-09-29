% l/vlPAG Pre Saccharine Consumption
% 0/15/2024


cd('/Users/nitishpatel/Library/CloudStorage/Box-Box/');

ch_on_off = readmatrix('PAG_VGAT_Feeding.xlsx','Sheet','Fig3_RevalueChow_DevalueSacchar','Range','M76:M81');
yfp_on_off = readmatrix('PAG_VGAT_Feeding.xlsx','Sheet','Fig3_RevalueChow_DevalueSacchar','Range','M91:M96');


p = ranksum(ch_on_off, yfp_on_off);

figure;
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
yticks([-25 150])
ylim([-25 150])
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

ylabel({'\Delta saccharine weight (mg)'; '(ON-OFF)'}, 'FontSize', 24);
text(0.24, 154, aster, 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
text(0.24, 160, '___', 'HorizontalAlignment', 'center', 'FontSize', 40);
