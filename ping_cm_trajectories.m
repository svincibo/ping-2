clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

measure = 'myelinmap'; % 'fa', 'volume', 'gmd', 'md'

hemisphere = 'both'; %left, right , both

save_figures = 'yes';
alphastat = 0.66; % to return 1 SD, for 95% CI use .05

color_adults = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
color_children = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed.
    %     outlier = [108 126 318];%
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = {''};
    
    
else
    
    outlier = {''};
    
end

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

if strcmp(measure, 'fa')
    ylimlo = 0.12; ylimhi = 0.25;
elseif strcmp(measure, 'volume')
    ylimlo = 0.001; ylimhi = .005;
    %             ylimlo = 0.0000006; ylimhi = .0000025;
elseif strcmp(measure, 'snr')
    ylimlo = 5; ylimhi = 32;
elseif strcmp(measure, 'gmd')
    ylimlo = 0.8; ylimhi = .95;
    elseif strcmp(measure, 'myelinmap')
    ylimlo = 1.5; ylimhi = 3.0;
elseif strcmp(measure, 'thickness')
    ylimlo = 1.5; ylimhi = 3.5;
end

%% WHITE MATTER MEASURES

% Read in data.
d = readtable(fullfile(rootDir, 'supportFiles', ['ping_cm_data_forSPSS_' measure '.csv']));

% Get index for outliers to be removed.
idx_keep = find(~ismember(d.subID, outlier));

% Remove outliers.
d = d(idx_keep, :);

scatter(d.age_mri, d.occipital, '.r'); hold on;
scatter(d.age_mri, d.ventral, '.g')
scatter(d.age_mri, d.parietal, '.b')
scatter(d.age_mri, d.frontal, '.k')


% Tell Matlab that sex and age group are categorical variables.
d.sex = categorical(d.sex);
d.age_mri = categorical(d.age_mri);
d.age_mri_tag = categorical(d.age_mri_tag);
d.age_beh = categorical(d.age_beh);
d.scanner_mri = categorical(d.scanner_mri);

% Get residuals after removing sex effects.
modelspec_confound = 'occipital ~ sex';
mdlr = fitlm(d, modelspec_confound);
d.res = table2array(mdlr.Residuals(:, 1));

% Display effect of confound.
disp(modelspec_confound)
disp(mdlr.anova)

%% LOBES.

roinames = {'occipital', 'ventral', 'parietal', 'frontal'};

for r = 1:length(roinames)
    
        % Mean center continuous variables for modelling.
    m_tract_demeaned = double(m_tract - nanmean(m_tract));
    
% 1. Linear main effects model.
    modelspec = 'res ~ age_mri';
            disp([roinames{r} ', ' modelspec])
    if sum(remove) == 0
        
        % Fit regression model.
        mdlr1 = fitlm(d, modelspec);
        
    else
        
        % Fit regression model, excluding outliers.
        %         mdlr1 = fitlm(data, modelspec, 'Exclude', find(sum(data.subID == remove, 2)));
        mdlr1 = fitlm(d, modelspec);
        
    end
    
    % Get model stats.
    stat = mdlr1.anova; disp(stat);
    
    % Write to file.
%     fprintf(fid, '%s \t 1 %s \t\t %1.3f \t %3.3f \t %1.3f \t %1.3f \n', roinames{r}, modelspec, ...
%         mdlr1.Rsquared.Adjusted, mdlr1.ModelCriterion.AICc, round(stat.F(2), 3), round(stat.pValue(2), 3));
   
    % 2. Nonlinear main effects model.
        modelspec = 'res ~ cov_age^2';
            disp([roinames{r} ', ' modelspec])
        if sum(remove) == 0
            
            % Fit regression model.
            mdlr2 = fitlm(d, modelspec);
            
        else
            
            % Fit regression model, excluding outliers.
            mdlr2 = fitlm(d, modelspec, 'Exclude', find(sum(d.subID == remove, 2)));
            
        end
        
        % Get model stats.
        stat = mdlr2.anova; disp(stat);
        
        % Write to file.
%         fprintf(fid, '%s \t 2 %s \t\t %1.3f \t %3.3f \t %1.3f \t %1.3f \n', roinames{r}, modelspec, ...
%             mdlr2.Rsquared.Adjusted, mdlr2.ModelCriterion.AICc, round(stat.F(2), 3), round(stat.pValue(2), 3));
        
end




occipital_child = d.occipital(find(contains(d.age_mri_tag, 'child_upto8.5')));
ventral_child = d.ventral(find(contains(d.age_mri_tag, 'child_upto8.5')));
parietal_child = d.parietal(find(contains(d.age_mri_tag, 'child_upto8.5')));
frontal_child = d.frontal(find(contains(d.age_mri_tag, 'child_upto8.5')));

occipital_adult = d.occipital(find(contains(d.age_mri_tag, 'adult')));
ventral_adult = d.ventral(find(contains(d.age_mri_tag, 'adult')));
parietal_adult = d.parietal(find(contains(d.age_mri_tag, 'adult')));
frontal_adult = d.frontal(find(contains(d.age_mri_tag, 'adult')));

%% Plot means including individual data.

figure(1)
hold on;

coloralpha = .1;

% Means (do this first for legend and then second to keep it on top layer).
child_mean = [nanmean(occipital_child) nanmean(ventral_child) nanmean(parietal_child) nanmean(frontal_child)];
adult_mean = [nanmean(occipital_adult) nanmean(ventral_adult) nanmean(parietal_adult) nanmean(frontal_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

% Individual data points for children.
scatter(repmat(1, size(occipital_child)), occipital_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(ventral_child)), ventral_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(parietal_child)), parietal_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(frontal_child)), frontal_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Individual data points for adults.
scatter(repmat(1, size(occipital_adult)), occipital_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(ventral_adult)), ventral_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(parietal_adult)), parietal_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(frontal_adult)), frontal_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Means (second time to put it on top layer).
child_mean = [nanmean(occipital_child) nanmean(ventral_child) nanmean(parietal_child) nanmean(frontal_child)];
adult_mean = [nanmean(occipital_adult) nanmean(ventral_adult) nanmean(parietal_adult) nanmean(frontal_adult)];

child_sd = [std(occipital_child, 'omitnan') std(ventral_child, 'omitnan') std(parietal_child, 'omitnan') std(frontal_child, 'omitnan')];
adult_sd = [std(occipital_adult, 'omitnan') std(ventral_adult, 'omitnan') std(parietal_adult, 'omitnan') std(frontal_adult, 'omitnan')];

% child_ci = 1.96*[std(occipital_child)/sqrt(vof_n_child) std(ventral_child)/sqrt(hv_n_child) std(parietal_child)/sqrt(pv_n_child) std(frontal_child)/sqrt(hd_n_child) std(fat_child)/sqrt(hd_n_child)];
% adult_ci = 1.96*[std(occipital_adult)/sqrt(vof_n_adult) std(ventral_adult)/sqrt(hv_n_adult) std(parietal_adult)/sqrt(pv_n_adult) std(frontal_adult)/sqrt(hd_n_adult) std(fat_adult)/sqrt(hd_n_adult)];

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
xlabels = {'Occipital', 'Ventral', 'Parietal', 'Frontal'};
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
if strcmp(measure, 'volume')
    yax.TickLabels = {num2str(ylimlo, '%2.1e'), num2str((ylimlo+ylimhi)/2, '%2.1e'), num2str(ylimhi, '%2.1e')};
else
    yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
end
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

% legend({'Children', 'Adults'}, 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', fontsize);
% legend('boxoff');

if strcmp(measure, 'fa')
    a.YLabel.String = 'Fractional Anisotropy (FA)';
elseif strcmp(measure, 'volume')
    a.YLabel.String = 'Gray Matter Volume (GMV) Proportion';
elseif strcmp(measure, 'gmd')
    a.YLabel.String = 'Gray Matter Density (GMD)';
elseif strcmp(measure, 'snr')
    a.YLabel.String = 'tSNR';
elseif strcmp(measure, 'myelinmap')
    a.YLabel.String = 't1/t2 ratio';
elseif strcmp(measure, 'thickness')
    a.YLabel.String = 'Thickness';
end

a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_cm_trajectories_' measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_cm_trajectories_' measure '_' hemisphere]), '-depsc')
    
end

hold off;



