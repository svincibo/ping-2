% Rearranges data in ping_data_fa.mat (output of
% ping_datacat_mri_beh_forSPSS.m) into a format that SPSS likes.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

% beh_measure = 'age'; %age, lit, vm, fm
wm_measure_here = {'fa'};
streamline_min = 100;

% % Load wm data from ping_qa_tractstats.m.
load(fullfile(rootDir, 'supportFiles', 'ping_data_streamlinecount.mat'));
streamlinecounts_subID = streamlinecounts(:, 1);
streamlinecounts = streamlinecounts(:, 4:end);
streamline_idx = table2array(streamlinecounts) >= streamline_min;

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data.
    load(fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure_here{w} '.mat']))
    d = data_all;
    clear data_all
    
    % Convert into array and header for ease.
    %     data_all_in = table2array(data_tbl);
    %     data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Select tois from d.
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices tracts of interest: col.
        t_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2') || strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftSLF3') || strcmp(d.Properties.VariableNames{k}, 'rightSLF3') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftAslant') || strcmp(d.Properties.VariableNames{k}, 'rightAslant') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftILF') || strcmp(d.Properties.VariableNames{k}, 'rightILF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftIFOF') || strcmp(d.Properties.VariableNames{k}, 'rightIFOF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftTPC') || strcmp(d.Properties.VariableNames{k}, 'rightTPC') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftpArc') || strcmp(d.Properties.VariableNames{k}, 'rightpArc') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftVOF') || strcmp(d.Properties.VariableNames{k}, 'rightVOF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftMDLFang') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFang');
        
    end
    t_temp = d(:, t_idx);
    
    % Change cells with low streamline counts to NaN.
    for k = 1:length(t_temp.Properties.VariableNames)
        
        % Reorganize streamline_idx into lostream_idx so that the columns are in the same order as toi_idx.
        col_idx = find(contains(streamlinecounts.Properties.VariableNames, t_temp.Properties.VariableNames{k}));
        
        if ~isempty(col_idx)
            
            t(:, k) = streamline_idx(:, col_idx).*table2array(t_temp(:, k));
            
        end
        
    end
    t(t == 0) = NaN;
    
    % Return to table format.
    toi = array2table(t, 'VariableNames', d.Properties.VariableNames(t_idx));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(table2array(toi(find(d.gp_age == 1), :)));
    m_wtwg3 = nanmean(table2array(toi(find(d.gp_age == 3), :)));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(table2array(toi(find(d.gp_age == 1), :)));
    std_wtwg3 = nanstd(table2array(toi(find(d.gp_age == 3), :)));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3;
    r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(table2array(toi) > r_max_wtwg1))
        disp(['Replaced ' num2str(numel(find(table2array(toi) > r_max_wtwg1))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.'])
        [idx, idy] = find(table2array(toi) > r_max_wtwg1);
        if ~isempty(idx)
            toi(idx, idy) = {NaN};
        end
    else
        disp('No data points were above35 standard deviations of the within-tract, within-group mean.')
    end
    %Min
    if ~isempty(find(table2array(toi) < r_min_wtwg1))
        disp(['Replaced ' num2str(numel(find(table2array(toi) < r_min_wtwg1))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.'])
        [idx, idy] = find(table2array(toi) < r_min_wtwg1);
        if ~isempty(idx)
            toi(idx, idy) = {NaN};
        end
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.')
    end
    
    %     % Make hemisphere averaged column for each tract.
    %     toi_havg = cat(2, nanmean(toi(:, [1 11]), 2), nanmean(toi(:, [2 12]), 2), nanmean(toi(:, [3 13]), 2), nanmean(toi(:, [4 14]), 2), nanmean(toi(:, [5 15]), 2), ...
    %         nanmean(toi(:, [6 16]), 2), nanmean(toi(:, [7 17]), 2), nanmean(toi(:, [8 18]), 2), nanmean(toi(:, [9 19]), 2), nanmean(toi(:, [10 20]), 2));
    %     toi_havg_header = erase(data_all_in_header(1:10), 'left');
    
    % Make tract group
    
    
    % Output csv file for ANOVA in SPSS. (Matlab doesn't handle Mixed Model
    % ANOVAs well when the between-group variable is correlated with subID
    % (e.g., when between-group variable is something like age groups).
    t_out = [table(d.subID, 'VariableNames', {'subID'}) table(d.group_str, 'VariableNames', {'group_str'}) ...
        table(d.cov_age, 'VariableNames', {'cov_age'}) table(d.cov_sex, 'VariableNames', {'cov_sex'}) ...
        table(d.read_beh, 'VariableNames', {'read_beh'}) table(d.scanner, 'VariableNames', {'scanner'}) ...
        table(d.gp_age, 'VariableNames', {'gp_age'}) table(d.gp_scanner, 'VariableNames', {'gp_scanner'}) toi];
    
    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure_here{w} '_forSPSS.csv']));
    fid = fopen(fullfile(rootDir, 'supportFiles', ['ping_data_' wm_measure_here{w} '_forSPSS.csv']));
    fclose(fid);
    
%     %     % Output z-scored file for SPSS.
%     d = d(d.gp_age~=3, :);
%     beh_z = (table2array(d(:, [5 7 9 11:19])) - nanmean(table2array(d(:, [5 7 9 11:19])), 1))./nanstd(table2array(d(:, [5 7 9 11:19])), [], 1);
%     toi = toi(d.gp_age ~=3, :);
%     toi_z = (table2array(toi) - nanmean(table2array(toi), 1))./nanstd(table2array(toi), [], 1);
% %     toi_havg_z = (nanmean(toi_havg, 1) - toi_havg)./nanstd(toi_havg, [], 1);
%     
%     %     temp = nanmean(toi(:, hv == 1), 2);
%     %     toi_meanhd_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
%     %
%     %     temp = nanmean(toi(:, hv == 2), 2);
%     %     toi_meanhv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
%     %
%     %     temp = nanmean(toi(:, hv == 3), 2);
%     %     toi_meanv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
%     
%     v = strcat('z_', toi.Properties.VariableNames(2:end));
%     t_out_z = array2table(cat(2, d.subID, d.gp_age, d.gp_lit, d.gp_vm, ...
%         d.gp_fm, d.cov_age, d.cov_sex, beh_z, toi_z(:, 2:end)), ...
%         'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', 'z_c_lit', 'z_c_vm','z_c_fm','z_BeeryVMI', 'z_BeeryVP', 'z_BeeryMC', 'z_gPegs_dom_forward', 'z_gPegs_nondom_forward', ...
%         'z_WJIV_LetterWordIdentification', 'z_WJIV_Spelling', 'z_WJIV_WordAttack', 'z_WJIV_SpellingOfSounds', v{:, :}});
%     writetable(t_out_z, fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell_z.csv']));
%     fid = fopen(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell_z.csv']));
%     fclose(fid)
%     
end


