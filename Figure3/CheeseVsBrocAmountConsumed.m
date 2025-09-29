cd('/Users/nitishpatel/Library/CloudStorage/Box-Box/'); 
data = [ readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','O8:O19'),... 
    readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','D8:D19'),... 
    readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','R8:R19'),... 
    readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','G8:G19')... 
    ]; 
    
data(5,:)  = []; 
data(10,:) = [];

data_onoff = [data(:,2)-data(:,1), data(:,4)-data(:,3)];
means = mean(data_onoff,1);
sems  = std(data_onoff,[],1) ./ sqrt(size(data_onoff,1));

% Colors
c1 = [0.2 0.5 0.8];   % blueish
c2 = [0.8 0.4 0.4];   % reddish

%% --------- Figure 1 (first column) ---------
figure; hold on;
bar(1, means(1), 0.3, 'FaceColor', c1, 'EdgeColor','none');
errorbar(1, means(1), sems(1), 'k', 'linestyle','none', 'LineWidth',1.2);
scatter( ones(size(data_onoff,1),1), data_onoff(:,1), ...
         40, 'ok', 'MarkerFaceColor', 'k');

xlim([0.5 1.5])
ylim([0 300])
yticks([0 300])
set(gca,'linewidth',2,'TickDir','out','xtick',[],'xticklabel',[], ...
    'XAxisLocation','origin','fontsize',30)
ylabel({'amount consumed (mg)';'(ON-OFF)'});
text(1, 300, "**", 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
box off; grid off;
%exportgraphics(gca,'cheese vs broc opto bars - cheese.eps','ContentType','vector')

%% --------- Figure 2 (second column) ---------
figure; hold on;
bar(1, means(2), 0.3, 'FaceColor', c2, 'EdgeColor','none');
errorbar(1, means(2), sems(2), 'k', 'linestyle','none', 'LineWidth',1.2);
scatter( ones(size(data_onoff,1),1), data_onoff(:,2), ...
         40, 'ok', 'MarkerFaceColor', 'k');

xlim([0.5 1.5])
ylim([-30 10])
yticks([-30 10])
set(gca,'linewidth',2,'TickDir','out','xtick',[],'xticklabel',[], ...
    'XAxisLocation','origin','fontsize',30)
text(1, 10, "*", 'HorizontalAlignment', 'center', 'FontSize', 30); % adds significant asterisk
ylabel({'amount consumed (mg)';'(ON-OFF)'});

box off; grid off;
%exportgraphics(gca,'cheese vs broc opto bars - broccoli.eps','ContentType','vector')


%%%%%%%%%%%%%%%%%%%%%%%%
% Pie charts
values = [64.73644221 35.26355779];
labels = {'Cheese','Broccoli'};

figure;
p = pie(values);

% Attach labels with percentages
percentLabels = strcat(labels, {' '}, compose('%.1f%%', 100*values/sum(values)));
legend(percentLabels, 'Location','eastoutside');
%exportgraphics(gca,'cheese vs broc laserOFF.eps','ContentType','vector')

values = [87.43543581 12.56456419];
labels = {'Cheese','Broccoli'};

figure;
p = pie(values);

% Attach labels with percentages
percentLabels = strcat(labels, {' '}, compose('%.1f%%', 100*values/sum(values)));
legend(percentLabels, 'Location','eastoutside');
%exportgraphics(gca,'cheese vs broc laserON.eps','ContentType','vector')


%%%%%%%%%%%%%%%%%%%%%%%
% Paired permutation test

% Paired permutation test for cheese preference (ON vs OFF)
cd('/Users/nitishpatel/Library/CloudStorage/Box-Box/'); 
% Example data (replace with your own)
OFF_cheese = readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','S8:S19');
ON_cheese  = readmatrix('PAG_VGAT_Feeding','Sheet','BrocVsCheese','Range','H8:H19');

OFF_cheese(5,:) = [];
OFF_cheese(10,:) = [];
ON_cheese(5,:) = [];
ON_cheese(10,:) = [];

OFF_broc = 1 - OFF_cheese;
ON_broc  = 1 - ON_cheese;

nIter = 10000;

% Observed differences
diffs_cheese = ON_cheese - OFF_cheese;
diffs_broc   = ON_broc  - OFF_broc;

obs_cheese = mean(diffs_cheese);
obs_broc   = mean(diffs_broc);

% Joint statistic: cheese increase AND broccoli decrease
obs_joint = min(obs_cheese, -obs_broc);

% Null distribution
null_joint = zeros(1,nIter);
for i = 1:nIter
    flipSigns = randi([0 1], size(diffs_cheese))*2 - 1;
    null_cheese = mean(diffs_cheese .* flipSigns);
    null_broc   = mean(diffs_broc   .* flipSigns);
    null_joint(i) = min(null_cheese, -null_broc);
end

% One-tailed test (we only care about cheese↑ and broccoli↓)
p_joint = mean(null_joint >= obs_joint);

% Results
fprintf('Observed cheese diff = %.3f\n', obs_cheese);
fprintf('Observed broc diff   = %.3f\n', obs_broc);
fprintf('Joint stat           = %.3f\n', obs_joint);
fprintf('Permutation p-value (cheese↑ & broccoli↓) = %.4f\n', p_joint);

% Plot null distribution
figure;
histogram(null_joint,'Normalization','pdf','FaceColor',[0.5 0.5 0.5]);
hold on;
yl = ylim;
plot([obs_joint obs_joint], yl, 'r-', 'LineWidth', 2);
xlabel('Joint statistic (cheese↑ & broccoli↓)');
ylabel('Probability density');
title(sprintf('Joint permutation test (p = %.4f)', p_joint));