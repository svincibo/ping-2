clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

wm_measure = 'fa';

hemisphere = 'both'; %left, right , both

save_figures = 'yes';
alphastat = 0.66; % to return 1 SD, for 95% CI use .05

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
    % P1039, P1246, P0787, P1409
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = {'P0133', 'P0201', 'P0222', 'P0688', 'P1012', 'P0815', 'P1039', 'P1246', 'P0787', 'P1409'};
      
else
    
    outlier = [];
    
end

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 50;
xtickvalues = [1 2 3 4 5];
xlim_lo = 0.5; xlim_hi = 5.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = 0.35; ylimhi = 0.60;

%% WHITE MATTER MEASURES

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
d = readtable(fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure '_forSPSS.csv']));

if strcmp(remove_outliers, 'yes')
    
    % Get index for outliers to be removed.
    idx_keep = find(~contains(d.subID, outlier));
    
    % Remove outliers.
    d = d(idx_keep, :);
    
end

% Get easy index for age group.
group = d.gp_age;

% Find children in the upto8.5 group that are not yet 5 or are above 8.
group(find((d.cov_age >= 4.5 & d.cov_age <= 8.5) & group == 1)) = 5; % just pick a number that means nothing so that the rest of the code won't recognize it as anything

% VOF
if strcmp(hemisphere, 'both')
    vof_child = nanmean(cat(2, d.leftVOF(group == 1), d.rightVOF(group == 1)), 2);
    vof_adult = nanmean(cat(2, d.leftVOF(group == 4), d.rightVOF(group == 4)), 2);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
elseif strcmp(hemisphere, 'left')
    vof_child = d.leftVOF(group == 1);
    vof_adult = d.leftVOF(group == 4);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
elseif strcmp(hemisphere, 'right')
    vof_child = d.rightVOF(group == 1);
    vof_adult = d.rightVOF(group == 4);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
end

% Horizontal, Ventral
if strcmp(hemisphere, 'both')
    hv_child = nanmean(cat(2, d.leftILF(group == 1), d.rightILF(group == 1), ...
        d.leftIFOF(group == 1), d.rightIFOF(group == 1)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 4), d.rightILF(group == 4), ...
        d.leftIFOF(group == 4), d.rightIFOF(group == 4)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(hemisphere, 'left')
    hv_child = nanmean(cat(2, d.leftILF(group == 1), d.leftIFOF(group == 1)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 4), d.leftIFOF(group == 4)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(hemisphere, 'right')
    hv_child = nanmean(cat(2, d.rightILF(group == 1), d.rightIFOF(group == 1)), 2);
    hv_adult = nanmean(cat(2, d.rightILF(group == 4), d.rightIFOF(group == 4)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
end

% % Horizontal, Ventral
% if strcmp(hemisphere, 'both')
%     hv_child = nanmean(cat(2, d.leftILF(group == 1), d.rightILF(group == 1)), 2);
%     hv_adult = nanmean(cat(2, d.leftILF(group == 4), d.rightILF(group == 4)), 2);
%     hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
% elseif strcmp(hemisphere, 'left')
%     hv_child = nanmean(d.leftILF(group == 1), 2);
%     hv_adult = nanmean(d.leftILF(group == 4), 2);
%     hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
% elseif strcmp(hemisphere, 'right')
%     hv_child = nanmean(d.rightILF(group == 1), 2);
%     hv_adult = nanmean(d.rightILF(group == 4), 2);
%     hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
% end

% Vertical, Posterior
if strcmp(hemisphere, 'both')
    pv_child = nanmean(cat(2, d.leftMDLFang(group == 1), d.rightMDLFang(group == 1), ...
        d.leftMDLFspl(group == 1), d.rightMDLFspl(group == 1), ...
        d.leftTPC(group == 1), d.rightTPC(group == 1), ...
        d.leftpArc(group == 1), d.rightpArc(group == 1)), 2);
    pv_adult = nanmean(cat(2, d.leftMDLFang(group == 4), d.rightMDLFang(group == 4), ...
        d.leftMDLFspl(group == 4), d.rightMDLFspl(group == 4), ...
        d.leftTPC(group == 4), d.rightTPC(group == 4), ...
        d.leftpArc(group == 4), d.rightpArc(group == 4)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
elseif strcmp(hemisphere, 'left')
    pv_child = nanmean(cat(2, d.leftMDLFang(group == 1), d.leftMDLFspl(group == 1), d.leftTPC(group == 1), d.rightTPC(group == 1), d.leftpArc(group == 1)), 2);
    pv_adult = nanmean(cat(2, d.leftMDLFang(group == 4), d.leftMDLFspl(group == 4), d.leftTPC(group == 4), d.leftpArc(group == 4)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
elseif strcmp(hemisphere, 'right')
    pv_child = nanmean(cat(2, d.rightMDLFang(group == 1), d.rightMDLFspl(group == 1), d.rightTPC(group == 1), d.rightTPC(group == 1), d.rightpArc(group == 1)), 2);
    pv_adult = nanmean(cat(2, d.rightMDLFang(group == 4), d.rightMDLFspl(group == 4), d.rightTPC(group == 4), d.rightpArc(group == 4)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
end

% Horizontal, Dorsal
if strcmp(hemisphere, 'both')
    hd_child = nanmean(cat(2, d.leftSLF1And2(group == 1), d.rightSLF1And2(group == 1), ...
        d.leftSLF3(group == 1), d.rightSLF3(group == 1)), 2);
    hd_adult = nanmean(cat(2, d.leftSLF1And2(group == 4), d.rightSLF1And2(group == 4), ...
        d.leftSLF3(group == 4), d.rightSLF3(group == 4)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
elseif strcmp(hemisphere, 'left')
    hd_child = nanmean(cat(2, d.leftSLF1And2(group == 1), d.leftSLF3(group == 1)), 2);
    hd_adult = nanmean(cat(2, d.leftSLF1And2(group == 4), d.leftSLF3(group == 4)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
elseif strcmp(hemisphere, 'right')
    hd_child = nanmean(cat(2, d.rightSLF1And2(group == 1), d.rightSLF3(group == 1)), 2);
    hd_adult = nanmean(cat(2, d.rightSLF1And2(group == 4), d.rightSLF3(group == 4)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
end

% FAT
if strcmp(hemisphere, 'both')
    fat_child = nanmean(cat(2, d.leftAslant(group == 1), d.rightAslant(group == 1)), 2);
    fat_adult = nanmean(cat(2, d.leftAslant(group == 4), d.rightAslant(group == 4)), 2);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
elseif strcmp(hemisphere, 'left')
    fat_child = d.leftAslant(group == 1);
    fat_adult = d.leftAslant(group == 4);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
elseif strcmp(hemisphere, 'right')
    fat_child = d.rightAslant(group == 1);
    fat_adult = d.rightAslant(group == 4);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
end

%% Plot means, but including individual data.
 
figure(1)
hold on;

coloralpha = .1;

% Means (do this first for legend and then second to keep it on top layer).
child_mean = [nanmean(vof_child) nanmean(hv_child) nanmean(pv_child) nanmean(hd_child) nanmean(fat_child)];
adult_mean = [nanmean(vof_adult) nanmean(hv_adult) nanmean(pv_adult) nanmean(hd_adult) nanmean(fat_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

% Individual data points for children.
scatter(repmat(1, size(vof_child)), vof_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(hv_child)), hv_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(pv_child)), pv_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(hd_child)), hd_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(5, size(fat_child)), fat_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Individual data points for adults.
scatter(repmat(1, size(vof_adult)), vof_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(hv_adult)), hv_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(pv_adult)), pv_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(hd_adult)), hd_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(5, size(fat_adult)), fat_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Means (second time to put it on top layer).
child_mean = [nanmean(vof_child) nanmean(hv_child) nanmean(pv_child) nanmean(hd_child) nanmean(fat_child)];
adult_mean = [nanmean(vof_adult) nanmean(hv_adult) nanmean(pv_adult) nanmean(hd_adult) nanmean(fat_adult)];

child_sd = [std(vof_child, 'omitnan') std(hv_child) std(pv_child) std(hd_child) std(fat_child, 'omitnan')];
adult_sd = [std(vof_adult, 'omitnan') std(hv_adult) std(pv_adult) std(hd_adult) std(fat_adult, 'omitnan')];

% child_ci = 1.96*[std(vof_child)/sqrt(vof_n_child) std(hv_child)/sqrt(hv_n_child) std(pv_child)/sqrt(pv_n_child) std(hd_child)/sqrt(hd_n_child) std(fat_child)/sqrt(hd_n_child)];
% adult_ci = 1.96*[std(vof_adult)/sqrt(vof_n_adult) std(hv_adult)/sqrt(hv_n_adult) std(pv_adult)/sqrt(pv_n_adult) std(hd_adult)/sqrt(hd_n_adult) std(fat_adult)/sqrt(hd_n_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

errorbar(xval, adult_mean, adult_sd, 'Color', color_adults, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
errorbar(xval, child_mean, child_sd, 'Color', color_children, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
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

legend({'Children', 'Adults'}, 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', fontsize);
legend('boxoff');

a.YLabel.String = 'Fractional Anisotropy (FA)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_fig3_ind_' wm_measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_fig3_ind_' wm_measure '_' hemisphere]), '-depsc')
    
end

hold off;

% Plot bootstrapped difference between children and adults.
figure(2)
hold on;

% Bootstrap to get diff and error bars because unequal sample sizes.
[vof_diff, vof_ci] = bootstrap_diff_unequalsamplesizes(vof_child, vof_adult, alphastat);
[hv_diff, hv_ci] = bootstrap_diff_unequalsamplesizes(hv_child, hv_adult, alphastat);
[pv_diff, pv_ci] = bootstrap_diff_unequalsamplesizes(pv_child, pv_adult, alphastat);
[hd_diff, hd_ci] = bootstrap_diff_unequalsamplesizes(hd_child, hd_adult, alphastat);
[fat_diff, fat_ci] = bootstrap_diff_unequalsamplesizes(fat_child, fat_adult, alphastat);

diff = [vof_diff hv_diff pv_diff hd_diff fat_diff];
sd = [vof_ci hv_ci pv_ci hd_ci fat_ci]; % returns 1 SD when alphastat = .66

xval = linspace(1, length(diff), length(diff));

scatter(xval, diff, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
errorbar(xval, diff, sd, 'Color', 'k', 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
ylimlo = -0.02; ylimhi = 0.10;
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo 0 (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), '0', num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

plot([xlim_lo xlim_hi], [0 0], 'k:')

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'Modularity Index (MI)';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-depsc')
    
end

hold off;

%% =========================================================================================== %%

function [diff_out, ci_out] = bootstrap_diff_unequalsamplesizes(group1, group2, alpha)

% get sample sizes of each group
n1 = length(group1);
n2 = length(group2);

% get distrution of differences
for r = 1:10000

    diff_real(r) = group2(randi(n2)) - group1(randi(n1));
     
end

diff_out = nanmean(diff_real);
ci_temp = prctile(diff_real, [100*alpha/2, 100*(1-alpha/2)]); % treats NaNs as missing values and removes them
ci_out = (ci_temp(2) - ci_temp(1));

end





