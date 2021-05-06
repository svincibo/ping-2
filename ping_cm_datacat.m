

clear all; close all; clc
format shortG

measure = {'thickness'}; % 'volume', 'fa', 'gmd', 'md', 'myelinmap', 'thickness'
roi = 'lobe'; %'lobe';

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

% Read in lobe codes.
lobe = readtable(fullfile(rootDir, 'supportFiles', 'LUT_glasser.csv'), 'TreatAsEmpty', {'.', 'na'});

% Read in parcellation codes.
parcel = readtable(fullfile(rootDir, 'supportFiles', 'key_glasser.csv'), 'TreatAsEmpty', {'.', 'na'});

% Read in behavioral data.
beh_in = readtable(fullfile(rootDir, 'supportFiles', 'ping_combined_participantinformation.csv'), 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'no';
if strcmp(remove_outliers, 'yes')
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = {'P0133'}; %P0133 is just an example, it is not to be removed
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

% % Parse table into array and header to make things easier.
% beh_data_in = table2array(beh);
% beh_data_in_header = beh.Properties.VariableNames;

for w = 1:length(measure)
    
    wm_measure = measure{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    scount = 0;
    % Load in measures for this subject.
    for i = 1:size(grp_contents, 1)
        
        %         if strcmp(measure, 'fa') || strcmp(measure, 'md')
        %
        %             % Get contents of the directory where the tract measures for this subject are stored.
        %             if strcmp(roi, 'lobe')
        %
        %                 sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-cortex_mapping_stats.id*', 'parc_MEAN.csv'));
        %
        %             elseif strcmp(roi, 'tractendpoints')
        %
        %                 sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-cortex_mapping_stats.tag-tract_endpoints*', 'tracts_MEAN.csv'));
        %
        %             end
        %
        %         elseif strcmp(measure, 'volume')
        %
        %             % Get contents of the directory where the tract measures for this subject are stored.
        %             sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-SupraTentorial*', 'rois.csv'));
        %
        %         elseif strcmp(measure, 'gmd')
        %
        %             % Get contents of the directory where the tract measures for this subject are stored.
        %             sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-gray_matter_density*', 'parc_MEAN.csv'));
        %
        %         elseif strcmp(measure, 'myelinmap')
        
        if strcmp(measure, 'myelinmap') || strcmp(measure, 'thickness')
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-myelin_mapping*', 'parc_MEAN.csv'));
            sub_contents_quality = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-myelin_mapping*', 'parc_COUNT_NONZERO.csv'));
        end
        %         end
        
        % Remove the '.' and '..' files.
        sub_contents_measure = sub_contents_measure(arrayfun(@(x) x.name(1), sub_contents_measure) ~= '.');
        sub_contents_quality = sub_contents_quality(arrayfun(@(x) x.name(1), sub_contents_quality) ~= '.');
        
        if ~isempty(sub_contents_measure)
            
            scount = scount + 1;
            
            % Grab subID.
            sub{scount} = grp_contents(i).name(end-4:end);
            
            % Display current sub ID.
            disp(grp_contents(scount).name)
            
            % Read in microstructural data for this subject.
            data_temp = readtable(fullfile(sub_contents_measure.folder, sub_contents_measure.name));
            
            % Read in quality data for this subject.
            quality_temp = readtable(fullfile(sub_contents_quality.folder, sub_contents_quality.name));
            
            % Change all measurements with less than 100 non-zero vertices to NaN.
            lo_quality_idx = find(table2array(quality_temp(:, 4)) < 500);
            data_temp{lo_quality_idx, 4:8} = NaN;
            
            %             % Change all measurements with snr <10 to NaN.
            %             lo_snr_idx = find(data_temp.snr<5);
            %             data_temp{lo_snr_idx, 4:8} = NaN;
            
            %             % Deal with idiosyncracies of the volume csv.
            %             if strcmp(measure, 'volume')
            %
            %                 % Remove the last 14 because those are subcortical.
            %                 data_temp = data_temp(1:end-14, :);
            %
            %                 % Rename ROI_names.
            %                 data_temp.ROI_name = str2double(erase(data_temp.ROI_name, 'ROI_'));
            %                 for r = 1:size(data_temp, 1)
            %
            %                     name_temp(r) = parcel.ROI_name(find(parcel.ROI_number == data_temp.ROI_name(r)));
            %
            %                 end
            %                 name_temp = erase(name_temp, {'lh.', 'rh.', '.label'});
            %                 data_temp.ROI_name = name_temp';
            %
            %             end
            
            % If data_out exists, append; if not, create.
            if i == 1
                
                % Create data_out array.
                data_out = data_temp;
                
            else
                
                % Concatenate this array with the previous subject's array.
                data_out = cat(1, data_out, data_temp);
                
            end
            
            clear data_temp
            
        end
        
        %         clear sub_contents_tractprofiles
        
    end % end i
    
    %     if strcmp(roi, 'tractendpoints')
    %
    %         % Select only ROIs from tracts that I care about.
    %         for t = 1:length(data_out.structureID)
    %
    %             idx(t) = contains(data_out.structureID{t}, 'leftSLF1And2') || contains(data_out.structureID{t}, 'rightSLF1And2') ...
    %                 || contains(data_out.structureID{t}, 'leftSLF3') || contains(data_out.structureID{t}, 'rightSLF3') ...
    %                 || contains(data_out.structureID{t}, 'leftAslant') || contains(data_out.structureID{t}, 'rightAslant') ...
    %                 || contains(data_out.structureID{t}, 'leftILF') || contains(data_out.structureID{t}, 'rightILF') ...
    %                 || contains(data_out.structureID{t}, 'leftIFOF') || contains(data_out.structureID{t}, 'rightIFOF') ...
    %                 || contains(data_out.structureID{t}, 'leftTPC') || contains(data_out.structureID{t}, 'rightTPC') ...
    %                 || contains(data_out.structureID{t}, 'leftpArc') || contains(data_out.structureID{t}, 'rightpArc') ...
    %                 || contains(data_out.structureID{t}, 'leftMDLFspl') || contains(data_out.structureID{t}, 'rightMDLFspl') ...
    %                 || contains(data_out.structureID{t}, 'leftVOF') || contains(data_out.structureID{t}, 'rightVOF') ...
    %                 || contains(data_out.structureID{t}, 'leftMDLFang') || contains(data_out.structureID{t}, 'rightMDLFang');
    %
    %         end
    %
    %         data_out = data_out(idx, :);
    %
    %     end
    
    %       % Get a list of unique cortex names.
    %     if strcmp(measure, 'fa') || strcmp(measure, 'md') || strcmp(measure, 'gmd')
    %         list_roi = unique(data_out.structureID);
    %     elseif strcmp(measure, 'volume')
    %         list_roi = unique(data_out.ROI_name);
    %     elseif strcmp(measure, 'myelinmap')
    list_roi = unique(data_out.structureID);
    %     end
    
    for y = 1:size(data_out, 1)
        
        name_temp = erase(data_out.structureID{y}, {'lh.', 'rh.', '_gaussian_1mm', '_FiberEndpoint'});
        data_out.structureID{y} = name_temp;
        
        if y == 1
            disp('Converting first 1/4 of ROI names ...');
        elseif y == ceil(0.25*size(data_out, 1))
            disp('Finished first 1/4. Converting the second 1/4 of ROI names ...');
        elseif y == ceil(0.5*size(data_out, 1))
            disp('Finished second 1/4. Converting the third 1/4 of ROI names ...');
        elseif y == ceil(0.75*size(data_out, 1))
            disp('Finished third 1/4. Converting the final 1/4 of ROI names ...');
        end
        
    end
    clear list_roi
    list_roi = unique(data_out.structureID);
    
    list_sub = unique(data_out.subjectID);
    
    % Get measurements for each ROI (reorganizing so that each column is an roi).
    for r = 1:size(list_roi, 1)
        
        for s = 1:length(list_sub)
            
            % Find index of this roi for this subject.
            idx_rs = find(strcmp(data_out.structureID, list_roi{r}) & strcmp(data_out.subjectID, list_sub{s}));
            
            if s == 1
                disp(['Starting ROI organization for ' data_out.structureID(idx_rs) '...']);
            end
            
            % Account for missing rois.
            if idx_rs
                
                if strcmp(measure, 'myelinmap')
                    measures(s, r) = data_out.myelinmap(idx_rs);
                elseif strcmp(measure, 'thickness')
                    measures(s, r) = data_out.thickness(idx_rs);
                end
                
                if r == 1
                    
                    disp(data_out.structureID(idx_rs));
                
                end
                
            else
                
                measures(s, r) = NaN;
                
            end
            
            clear idx_rs
            
        end
        
    end
    
    % Convert all zeros to NaN;
    measures(measures==0) = NaN;
    
    % Append subID.
    temp = table(transpose(sub), 'VariableNames', {'subID'});
    temp2 = array2table(measures, 'VariableNames', transpose(list_roi));
    wm = [temp temp2];
    clear temp temp2
    
    % Create grouping and behavioral vectors.
    beh = table(beh_in.subID, beh_in.sex, beh_in.age_mri, beh_in.age_mri_tag, beh_in.scanner_mri, beh_in.age_beh, beh_in.read_beh, ...
        'VariableNames', {'subID', 'sex', 'age_mri', 'age_mri_tag', 'scanner_mri', 'age_beh', 'read_beh'});
    
    % Join wm and beh tables based on subID, including only subIDs that are contained in both.
    data = join(wm, beh, 'Keys', 'subID');
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        if contains(data.subID, outlier)
            
            % Get index for outliers to be removed.
            idx_outlier = find(contains(data.subID, outlier));
            
            % Remove outliers.
            data = data(~idx_outlier, :);
            
        else
            
            data = data;
            
        end
        
    end
    
    %     % Create the output table.
    %     data_tbl = array2table(cat(2, data));
    
    % Save all variables.
    save(fullfile(rootDir, 'supportFiles', ['ping_cm_data_' wm_measure '.mat']), 'data')
    
    % Save as a CSV file.
    writetable(data, fullfile(rootDir, 'supportFiles', ['ping_cm_data_' wm_measure '.csv']))
    
    % Reset for next loop.
    %     clearvars -except w rootDir beh_data_in_tbl beh_data_in_header beh blprojectid remove_outliers w_measures outlier
    
end
