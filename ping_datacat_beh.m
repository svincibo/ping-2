% Read in csv files provided by PING for the mri sessions and the
% behavioral sessions (particularly the reading test sessions) and merge
% them into one csv containing only subjects that contain diffusion data
% and appear in both of the original csv and, also, including only
% variables of interest for the current study.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/Seagate/ping';

% Get bl project foldername.
blprojectid = 'proj-5fd37e5cdecf133f0daca92a';

% Read in demographic data from mri session and keep only participants for whom
% we have diffusion data. Note: This will also remove the first row that is the
% definition of the column variable.
mri = readtable(fullfile(rootDir, 'supportFiles', 'ping_mri_participantinformation.csv'), 'TreatAsEmpty', {'.', 'na'});
idx = find(contains(mri.image_description, 'diffusion'));
mri = mri(idx, :);

% Read in data from behavioral session (particularly the Reading Task data)
% and remove the first row that is the definition of the column variable.
beh = readtable(fullfile(rootDir, 'supportFiles/ping_beh_reading_participantinformation.csv'), 'TreatAsEmpty', {'.', 'na'});
beh = beh(2:end, :);

% Step through subjects in demo, one at a time, and construct a new table
% containing variables of interest.
s = 0;
for count = 1:size(mri, 1)
    
    % Do not include subjects that are to be removed for whatever reason (e.g., unexplainable redundancies).
    if ~strcmp(mri.src_subject_id{count}, 'P0008') && ~strcmp(mri.src_subject_id{count}, 'PING_PG_P1460')
        
        s = s + 1;
        
        % Find location of this suject in demo and orrt.
        mri_idx = find(contains(mri.src_subject_id, mri.src_subject_id{count}));
        beh_idx = find(contains(beh.src_subject_id, mri.src_subject_id{count}));
        
        % Gather the variables of interest.
        subID{s} = mri.src_subject_id{mri_idx};
        age_mri(s) = mri.interview_age(mri_idx);
        sex_mri{s} = mri.sex{mri_idx};
        scanner_mri{s} = mri.scanner_type_pd{mri_idx};
        
        age_beh(s) = beh.interview_age(beh_idx);
        sex_beh{s} = beh.sex{beh_idx};
        read_beh(s) = beh.tbx_reading_score(beh_idx);
        
        % Make needed tags.
        if age_mri(s)/12 <= 8.5
            
            age_mri_tag{s} = 'child_upto8.5';
            
        elseif age_mri(s)/12 > 8.5 && age_mri(s)/12 <= 13.5
            
            age_mri_tag{s} = 'child_8.6-13.5';
            
        elseif age_mri(s) > 13.5 && age_mri(s)/12 <= 18.5
            
            age_mri_tag{s} = 'child_13.6-18.5';
            
        elseif age_mri(s)/12 > 18.5
            
            age_mri_tag{s} = 'adult';
            
        else
            
            age_mri_tag{s} = 'NA';
            
        end
        
    end % end if ~strcmp
    
end % end s

% Check for inconsistencies in sex and output table accordingly.
if isequal(sex_mri, sex_beh)
    
    disp('No inconsistencies between mri and beh participant information regarding the variable sex.');
    
    % Combine all variables of interest into a new table.
    t_out = table(subID', sex_mri', age_mri', age_mri_tag', scanner_mri', age_beh', read_beh', ...
        'VariableNames', {'subID', 'sex', 'age_mri', 'age_mri_tag', 'scanner_mri', 'age_beh', 'read_beh'});
    
    % Write table.
    writetable(t_out, fullfile(rootDir, 'supportFiles/ping_combined_participantinformation.csv'));
    
    % Write out mat file.
    save(fullfile(rootDir, 'supportFiles/ping_combined_participantinformation.mat'), 't_out');
    
else
    
    disp('Please check output csv for inconsistencies between mri and beh participant information regarding the variable sex.');
    
    % Combine all variables of interest into a new table.
    t_out = table(subID', sex_mri', sex_beh', 'VariableNames', {'subID', 'sex_mri', 'sex_beh'});
    
    % Write table.
    writetable(t_out, fullfile(rootDir, 'supportFiles/ping_checkinconsistenciesbetweenmriandbehsex.csv'));
    
end


