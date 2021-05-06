clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

wm_measure = 'fa';
hemisphere = 'both';
save_figures = 'yes';

% Specify colors.
ifof = [142 198 255]/255; % light blue
ilf = [0 127 255]/255; % dark blue

slf12 = [237 177 32]/255; % burnt yellow
slf3 = [240 221 165]/255; % light burnt yellow

parc = [64 224 208]/255; % turquoise
tpc =  [27 102 87]/255; % dark turquoise
mdlfspl = [42, 102, 0]/255; % green
mdlfang =  [207 255 226]/255; % sea foam

vof = [147 112 219]/255; % medium purple

fat = [240 128 128]/255; % light coral

color_adults = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
color_children = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed.
    % Identify outliers to be removed.
    % P0222, P0201, P0133, snr is below 20
    % P0688, streamline count too low in all right VHP tracts to get an average
    % P1012, P0815, streamline count too low in all leftPVP tracts to get an average
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = {'P0133', 'P0201', 'P0222', 'P0688', 'P1012', 'P0815'};
    
else
    
    outlier = {''};
    
end

% Set up plot and measure-specific details.
coloralpha = 0.2;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 50;
xlim_lo = 0; xlim_hi = 20;
xtickvalues = xlim_lo:5:xlim_hi;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = 0.3; ylimhi = 0.6;

%% WHITE MATTER MEASURES

d = readtable(fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure '_forSPSS.csv']));

if strcmp(remove_outliers, 'yes')
    
    % Get index for outliers to be removed.
    idx_keep = find(~contains(d.subID, outlier));
    
    % Remove outliers.
    d = d(idx_keep, :);
    
end

% Convert age from months to years.
d.cov_age = d.cov_age/12;

% Convert sex to variable.
females = find(strcmp(d.cov_sex, 'F')); d.cov_sex2(females) = 1;
males = find(strcmp(d.cov_sex, 'M')); d.cov_sex2(males) = 2;

% Make tract pathway averages.
if strcmp(hemisphere, 'both')
    d.vhp = nanmean([d.leftILF d.rightILF d.leftIFOF d.rightIFOF], 2);
    d.pvp = nanmean([d.leftMDLFang d.leftMDLFspl d.leftTPC d.leftpArc d.rightMDLFang d.rightMDLFspl d.rightTPC d.rightpArc], 2);
    d.dhp = nanmean([d.leftSLF1And2 d.leftSLF3 d.rightSLF1And2 d.rightSLF3], 2);
elseif strcmp(hemisphere, 'left')
    d.vhp = nanmean([d.leftILF d.leftIFOF], 2);
    d.pvp = nanmean([d.leftMDLFang d.leftMDLFspl d.leftTPC d.leftpArc], 2);
    d.dhp = nanmean([d.leftSLF1And2 d.leftSLF3], 2);
elseif strcmp(hemisphere, 'right')
    d.vhp = nanmean([d.rightILF d.rightIFOF], 2);
    d.pvp = nanmean([d.rightMDLFang d.rightMDLFspl d.rightTPC d.rightpArc], 2);
    d.dhp = nanmean([d.rightSLF1And2 d.rightSLF3], 2);
end

% Add box to highlight 4.5 - 8.5 years.
x = [4.5 8.5 8.5 4.5];
y = [0 0 1 1];
p2 = patch(x, y, [128 128 128]/255);
set(p2, 'facecolor', [128 128 128]/255, 'edgecolor', 'none', 'facealpha', .2); hold on;
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Scatter plots for tract categories
scatter(d.cov_age, d.vhp, ...
    'Marker', marker, 'SizeData', markersize, 'MarkerFaceColor', ilf, 'MarkerEdgeColor', ilf, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha*2); 
scatter(d.cov_age, d.pvp, ...
    'Marker', marker, 'SizeData', markersize, 'MarkerFaceColor', tpc, 'MarkerEdgeColor', tpc, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha*2);
scatter(d.cov_age, d.dhp, ...
    'Marker', marker, 'SizeData', markersize, 'MarkerFaceColor', slf12, 'MarkerEdgeColor', slf12, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha*2);

%% Get non-linear fit.

degp = 2;
x = d.cov_age;

% Tell Matlab that sex is categorical.
d.cov_sex2 = categorical(d.cov_sex2);

% VHP.
clr = ilf;
y = d.vhp;

% Specify the model.
g = fittype('a-b*exp(-c*x)');
f1 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
x1 = linspace(0, max(x), 1000);
plot(x1, f1(x1),'LineWidth', 3, 'LineStyle', '-', 'Color', clr);
% hi = f1 + std(f1); lo = f1 - std(f1); x2 = (1:size(f1, 2))';
% hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], clr);
% set(hp3, 'facecolor', clr, 'edgecolor', 'none', 'facealpha', .2);

% PVP.
clr = tpc;
y = d.pvp;

% % Plot fit.
g = fittype('a-b*exp(-c*x)');
f1 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
x1 = linspace(0, max(x), 1000);
plot(x1, f1(x1),'LineWidth', 3, 'LineStyle', '-', 'Color', clr);
% plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', clr)
% hi = f1 + std(f1); lo = f1 - std(f1); x2 = (1:size(f1, 2))';
% hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], clr);
% set(hp3, 'facecolor', clr, 'edgecolor', 'none', 'facealpha', .2);

% DHP.
clr = slf12;
y = d.dhp;

% Plot fit.
g = fittype('a-b*exp(-c*x)');
f1 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
x1 = linspace(0, max(x), 1000);
plot(x1, f1(x1),'LineWidth', 3, 'LineStyle', '-', 'Color', clr);
% plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', clr)
% hi = f1 + std(f1); lo = f1 - std(f1); x2 = (1:size(f1, 2))';
% hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], clr);
% set(hp3, 'facecolor', clr, 'edgecolor', 'none', 'facealpha', .2);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
% xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'VHP', 'PVP', 'DHP'}, 'Location', 'southeast', 'Orientation', 'vertical', 'FontSize', fontsize);
legend('boxoff');

a.YLabel.String = 'Fractional Anisotropy (FA)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
a.XLabel.String = 'Age (years)';
a.XLabel.FontSize = fontsize;
a.XLabel.FontAngle = 'normal';

pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_trajectories_' wm_measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_trajectories_' wm_measure '_' hemisphere]), '-depsc')
    
end

hold off;
