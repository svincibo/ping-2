% Reads in behavioral and mri measures and outputs a combined csv file and mat file.

clear all; close all; clc
format shortG

w_measures = {'fa'};%, 'md'};

fontname = 'Arial';
fontsizex = 16; fontsizey = 12;
fontangle = 'italic';
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
linewidth = .15;
alpha = .5;
save_figures = 'yes';

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

% Get bl project foldername.
blprojectid = 'proj-60708cf9c7f80a684995e0b1';

% Read in behavioral data.
beh = readtable(fullfile(rootDir, 'supportFiles', 'ping_combined_participantinformation.csv'), 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'no';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed.
    %     outlier = [108 126 318];%
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

% Parse table into array and header to make things easier.
% beh_data_in = table2array(beh);
beh_data_in_header = beh.Properties.VariableNames;

% Identify outliers to be removed.
% 128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
% 315 because strange WM in occipital lobe leading to no right or left VOF
% 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
% outlier = [128 315 318];

for w = 1:length(w_measures)
    
    wm_measure = w_measures{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    for i = 1:size(grp_contents, 1)
        
        % Grab subID.
        sub{i} = grp_contents(i).name(end-4:end);
        
        % Display current sub ID.
        disp(grp_contents(i).name)
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '/dt-neuro-tractprofile*/profiles/*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        for j = 1:size(sub_contents_tractprofiles)
            
%             % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
%             if i == 1 && j == 1
%                 
%                 tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
%                 
%             end
            
            % Read in data for this subject and this tract.
            data_temp = readtable([sub_contents_tractprofiles(j).folder filesep sub_contents_tractprofiles(j).name]);
            
            % Get middle 80%.
            start = size(data_temp, 1)*.1;
            stop = size(data_temp, 1)*.9;
            
            % Read in mean WM measure.
            if strcmp(wm_measure, 'ad')
                
                m(i, j) = nanmean(data_temp.ad_mean(start:stop));
                
            elseif strcmp(wm_measure, 'fa')
                
                m(i, j) = nanmean(data_temp.fa_mean(start:stop));
                
            elseif strcmp(wm_measure, 'md')
                
                m(i, j) = nanmean(data_temp.md_mean(start:stop));
                
            elseif strcmp(wm_measure, 'rd')
                
                m(i, j) = nanmean(data_temp.rd_mean(start:stop));
                
            end
            
            % Grab tract name for grouping variable.
            tract{i, j} = sub_contents_tractprofiles(j).name(1:end-13);
            
            clear data_temp
            
        end % end j
        
    end % end i
    
    % Find empty cells.
    t = find(cellfun(@isempty,tract));
    
    % Enter 'empty' in empty cells.
    tract(t) = {'empty'};
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get WM measurements for each tract (reorganizing so that each column is a tract).
    for k = 1:size(list_tract, 1)
        
        % Select the wm_measurements for this tract from each subject.
        temp = m.*contains(tract, list_tract{k});
        
        % Convert all zeros to NaN;
        temp(temp == 0) = NaN;
        
        % Get the mean of the wm_measure for this tract (take sum because don't want to include zero columns; only one value of interest per row).
        wm(:, k) = nansum(temp, 2);
        
        clear temp
        
    end % end k
    
    % Convert all zeros to NaN;
    wm(wm == 0) = NaN;
    
    % Remove 'empty' column from data and header and append subID.
    temp1 = table(transpose(sub), 'VariableNames', {'subID'});
    temp2 = array2table(wm(:, find(~all(isnan(wm), 1))));

    wm = horzcat(temp1, temp2);
    wm.Properties.VariableNames = [{'subID'}, transpose(list_tract(~strcmp(list_tract, 'empty')))];
    clear temp1 temp2
    
    % Create grouping and behavioral vectors.
    beh_out = table(beh.subID, beh.age_mri_tag, beh.age_mri, beh.sex, beh.read_beh, beh.scanner_mri);
    beh_out.Properties.VariableNames = {'subID', 'group_str', 'cov_age', 'cov_sex', 'read_beh', 'scanner'};
    
    % Append age grouping variable that is double.
    idx_child = contains(beh_out.group_str, 'child_upto8.5');
    idx_preteen = contains(beh_out.group_str, 'child_8.6-13.5');
    idx_adolescent = contains(beh_out.group_str, 'child_13.6-18.5');
    idx_adult = contains(beh_out.group_str, 'adult');
    beh_out.gp_Age = idx_child + 2*idx_preteen + 3*idx_adolescent + 4*idx_adult;
    
    % Append a scanner grouping variable that is double.
    idx_signa = contains(beh_out.scanner, 'SIGNA HDx');
    idx_trio = contains(beh_out.scanner, 'TrioTim');
    idx_discovery = contains(beh_out.scanner, 'DISCOVERY MR750');
    beh_out.gp_scanner = idx_signa + 2*idx_trio + 3*idx_discovery;
    
    % Determine which subIDs appear in both WM and BEH.
    sub_wm_beh = intersect(wm(:, find(strcmp(wm.Properties.VariableNames, 'subID'))), table(beh.subID, 'VariableNames', {'subID'}));
        
    % Get indices of subjects who appear in both WM and BEH.
    sub_idx_wm = ismember(wm(:, find(strcmp(wm.Properties.VariableNames, 'subID'))), sub_wm_beh);
    sub_idx_beh = ismember(table(beh.subID, 'VariableNames', {'subID'}), sub_wm_beh);

    % Select only subjects who appear in both WM and BEH.
    % Concatenate into one data array and one header array.
    % Remove redundant subID columns.
    data_all = [beh_out(sub_idx_beh, :) wm(sub_idx_wm, find(strcmp(wm.Properties.VariableNames, 'subID'))+1:end)];    
    data_all.Properties.VariableNames = [{'subID', 'group_str', 'cov_age', 'cov_sex', 'read_beh', 'scanner', 'gp_age', 'gp_scanner'}, ...
        wm.Properties.VariableNames{find(strcmp(wm.Properties.VariableNames, 'subID'))+1:end}];
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(data_all(:, find(strcmp(data_all.Properties.VariableNames, {'subID'}))), outlier);
        
        % Remove outliers.
        data_all = data_all(~idx_outlier, :);
        
    end
       
    % Save all variables.
    save(fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure '.mat']), 'data_all')
   
    % Write out table.
    writetable(data_all, fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure '.csv']));

    % Reset for next loop.
%     clearvars -except w rootDir beh_data_in_tbl beh_data_in_header beh blprojectid remove_outliers w_measures outlier
    
end
