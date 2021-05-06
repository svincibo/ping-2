% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

% beh_measure = 'age'; %age, lit, vm, fm
measure = {'thickness'};
roi = 'lobe'; %'tractendpoints';
streamline_min = 100;

% Read in lut for glasser atlas.
if strcmp(roi, 'lobe')
    lobe = readtable(fullfile(rootDir, 'supportFiles', 'LUT_glasser.csv'));
elseif strcmp(roi, 'tractendpoints')
    lobe = readtable(fullfile(rootDir, 'supportFiles', LUT_tractendpoints.csv'));
end

occipital = lobe.Name(lobe.Lobe == 1);
ventral = lobe.Name(lobe.Lobe == 2);
parietal = lobe.Name(lobe.Lobe == 3);
frontal = lobe.Name(lobe.Lobe == 4);

%% WHITE MATTER MEASURES
for w = 1:length(measure)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load(fullfile(rootDir, 'supportFiles', ['ping_cm_data_' measure{w} '.mat']))
    d = data;
    clear data
    
    % Convert into array and header for ease.
    %     data_all_in = table2array(data_tbl);
    %     data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Select tois from d.
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices tracts of interest: col.
        t_idx(k) = ismember(d.Properties.VariableNames{k}, occipital) || ismember(d.Properties.VariableNames{k}, ventral) ...
            || ismember(d.Properties.VariableNames{k}, parietal) ... 
            || ismember(d.Properties.VariableNames{k}, frontal) || strcmp(d.Properties.VariableNames{k}, 'icv');
%         || strcmp(d.Properties.VariableNames{k}, 'subID') 
        
    end
    roi = d(:, t_idx);
    
    % Return to table format.
    %     roi = array2table(t, 'VariableNames', d.Properties.VariableNames(t_idx));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(table2array(roi(find(contains(d.age_mri_tag, 'child_upto8.5')), :)));
    m_wtwg2 = nanmean(table2array(roi(find(contains(d.age_mri_tag, 'child_8.6-13.5')), :)));
    m_wtwg3 = nanmean(table2array(roi(find(contains(d.age_mri_tag, 'child_13.6-18.5')), :)));
    m_wtwg4 = nanmean(table2array(roi(find(contains(d.age_mri_tag, 'adult')), :)));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(table2array(roi(find(contains(d.age_mri_tag, 'child_upto8.5')), :)));
    std_wtwg2 = nanstd(table2array(roi(find(contains(d.age_mri_tag, 'child_8.6-13.5')), :)));
    std_wtwg3 = nanstd(table2array(roi(find(contains(d.age_mri_tag, 'child_13.6-18.5')), :)));
    std_wtwg4 = nanstd(table2array(roi(find(contains(d.age_mri_tag, 'adult')), :)));

    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1; r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2; r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3; r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    r_max_wtwg4 = m_wtwg4 + 3*std_wtwg4; r_min_wtwg4 = m_wtwg4 - 3*std_wtwg4;

    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, size(find(contains(d.age_mri_tag, 'child_upto8.5')))), ....
        repmat(r_max_wtwg2, size(find(contains(d.age_mri_tag, 'child_8.6-13.5')))), ...
        repmat(r_max_wtwg3, size(find(contains(d.age_mri_tag, 'child_13.6-18.5')))), ...
        repmat(r_max_wtwg4, size(find(contains(d.age_mri_tag, 'adult')))));
    r_min = cat(1, repmat(r_min_wtwg1, size(find(contains(d.age_mri_tag, 'child_upto8.5')))), ....
        repmat(r_min_wtwg2, size(find(contains(d.age_mri_tag, 'child_8.6-13.5')))), ...
        repmat(r_min_wtwg3, size(find(contains(d.age_mri_tag, 'child_13.6-18.5')))), ...
        repmat(r_min_wtwg4, size(find(contains(d.age_mri_tag, 'adult')))));
    
    % Display.
    disp([measure{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(table2array(roi) > r_max))
        disp(['Replaced ' num2str(numel(find(table2array(roi) > r_max))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.']);
        [idx, idy] = find(table2array(roi) > r_max);
        if ~isempty(idx)
            roi(idx, idy) = {NaN};
        end
    else
        disp('No data points were above 3 standard deviations of the within-tract, within-group mean.');
    end
    %Min
    if ~isempty(find(table2array(roi) < r_min))
        disp(['Replaced ' num2str(numel(find(table2array(roi) < r_min))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.']);
        [idx, idy] = find(table2array(roi) < r_min);
        if ~isempty(idx)
            roi(idx, idy) = {NaN};
        end
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.');
    end
    
    % Create lobe averages.
    for k = 1:length(roi.Properties.VariableNames)
        
        o_idx(k) = ismember(roi.Properties.VariableNames{k}, occipital);
        v_idx(k) = ismember(roi.Properties.VariableNames{k}, ventral);
        p_idx(k) = ismember(roi.Properties.VariableNames{k}, parietal);
        f_idx(k) = ismember(roi.Properties.VariableNames{k}, frontal);
        
    end
    
    % Get mean in each cortical region.
    odata = nanmean(table2array(roi(:, o_idx)), 2);
    vdata = nanmean(table2array(roi(:, v_idx)), 2);
    pdata = nanmean(table2array(roi(:, p_idx)), 2);
    fdata = nanmean(table2array(roi(:, f_idx)), 2);
    
%     % Creat lobe snr averages.
%     for k = 1:length(d.Properties.VariableNames)
%         
%         o_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(occipital, '_snr'));
%         v_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(ventral, '_snr'));
%         p_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(parietal, '_snr'));
%         f_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(frontal, '_snr'));
%         
%     end
%     
%     % Get mean in each cortical region.
%     osnrdata = nanmean(table2array(d(:, o_snr_idx)), 2);
%     vsnrdata = nanmean(table2array(d(:, v_snr_idx)), 2);
%     psnrdata = nanmean(table2array(d(:, p_snr_idx)), 2);
%     fsnrdata = nanmean(table2array(d(:, f_snr_idx)), 2);
%     
%             osnr = nanmean(table2array(d(:, o_snr_idx)), 2);
%         vsnr = nanmean(table2array(d(:, v_snr_idx)), 2);
%         psnr = nanmean(table2array(d(:, p_snr_idx)), 2);
%         fsnr = nanmean(table2array(d(:, f_snr_idx)), 2);
       
       
        % Output csv file for ANOVA in SPSS. (Matlab doesn't handle Mixed Model
        % ANOVAs well when the between-group variable is correlated with subID
        % (e.g., when between-group variable is something like age groups).
        temp = table(d.subID, d.sex, d.age_mri_tag, d.scanner_mri, 'VariableNames', {'subID', 'sex', 'age_mri_tag', 'scanner_mri'});
        temp2 = array2table(cat(2, d.age_mri, d.age_beh, d.read_beh, table2array(roi), odata, vdata, pdata, fdata), ...
            'VariableNames', {'age_mri', 'age_beh', 'read_beh', roi.Properties.VariableNames{:}, 'occipital', 'ventral', 'parietal', 'frontal'});
        t_out = [temp temp2];
        clear temp temp2

    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['ping_cm_data_forSPSS_' measure{w} '.csv']));
    fid = fopen(fullfile(rootDir, 'supportFiles', ['ping_cm_data_forSPSS_' measure{w} '.csv']));
    fclose(fid);
    
end


