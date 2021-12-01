function subj = roi_nifti_load_2(exp, sbj_id, session, roi_type, roi_name, varargin)
% makes subj structure and load nifti files to creat multi-voxel pattern
% -------------------------------------------------------------------------
% input:
%     exp       (string)    - experiment name (Lee2014)
%     sbj_id    (num)       - subject ID
%     session   (num)       - imaging session
%     roi_type  (string)    - type of ROI (AAL, AICHAmc, sphere10)
%     roi_name  (string)    - name of ROI 
%                 (lilPFC, rilPFC, F3TL, F3TR, V1L, V1R, HIPPOL, HIPPOR,
%                 VentStrL, VentStrR, )
% 
% output: subj (princeton-mvpa-toolbox-master)
% subj fields
%     mask      - ROI mask
%     pattern   - masked multi-voxel pattern
%     ...
% -------------------------------------------------------------------------
% <Name-value pairs>
% is_cluster
% preproc_prefix
% mask_num
% hemisphere
% save_dir


%% Setting & initialization

% default
options = struct('is_cluster', 1, ...
                'preproc_prefix', 'wra*', ...
                'mask_num', [], ...
                'hemisphere', [], ...
                'save_dir', []);
option_names = fieldnames(options);
if mod(length(varargin), 2) == 1
    error('(roi_nifti_load) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(roi_nifti_load) %s is not a recognized parameter name', pair{1})
    end
end

% root
if options.is_cluster
    root = '/home';
else
    root = '//143.248.30.94/bmlsamba';
end

% experiment setting
switch exp
    
    case 'Lee2014'
        % original subject IDs
        IDs = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24];
        
        % directories
        atlas_dir = [root '/ydsung/A_Research/princeton-mvpa-toolbox-master/ROI_mask/atlas'];
%         atlas_dir = 'Z:/JR/dimensionality/atlas';
        if options.is_cluster
            nifti_dir = [root sprintf('/ydsung/2014fmri/fmri_arbitration/od-arbitration-%03d', IDs(sbj_id))];
        else
            nifti_dir = ['D:/fmri_arbitration/' sprintf('od-arbitration-%03d', IDs(sbj_id))];
        end
        session_dir = sprintf('/func/run_%04d', session);
        
        
    case 'Kim2019'
        
        IDs = 1:24; IDs([9, 12, 16]) = [];
        
        % directories
%         atlas_dir = [root '/ydsung/A_Research/princeton-mvpa-toolbox-master/ROI_mask/atlas/' exp];
        % 2021-05-05
        atlas_dir = [root '/ydsung/A_Research/princeton-mvpa-toolbox-master/ROI_mask/atlas'];
        nifti_dir = [root sprintf('/ydsung/complexity/od-cog-%03d', IDs(sbj_id))];
        session_dir = sprintf('/func/run_%04d', session);
        
end

% subj structure initialization
subj = init_subj(exp, ['subject', num2str(sbj_id)]);



%% load ROI mask (.nii file)

switch roi_type
    
    case 'AAL'
        subj = load_spm_mask_ydsung(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/aal_2014epi_coord.nii'], ...
            'mask_num', options.mask_num, ...
            'hemisphere', options.hemisphere);
        
    case 'AAL3'
        
        switch roi_name
            case 'F3TL'
                mask_num = 9;
            case 'F3TR'
                mask_num = 10;
            case 'V1L'
                mask_num = 47;
            case 'V1R'
                mask_num = 48;                
        end
        
        if ~isempty(options.mask_num)
            mask_num = options.mask_num;
        end
        
        subj = load_spm_mask_ydsung(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/ROI_MNI_V7_2014epi_coord.nii'], ...
            'mask_num', mask_num, ...
            'hemisphere', options.hemisphere);
    
    case 'AICHAmc'
        
        switch roi_name
            case 'G_FIT'
                mask_num = 17;
        end
        
        if ~isempty(options.mask_num)
            mask_num = options.mask_num;
        end
        
        subj = load_spm_mask_ydsung(subj, [roi_type '_' roi_name], ...
            [atlas_dir '/AICHAmc_2014epi_coord.nii'], ...
            'mask_num', mask_num, ...
            'hemisphere', options.hemisphere);
        
end



%% load fMRI data (.nii file)

% file names
count = 0;
listing_s = dir([nifti_dir, session_dir, '/', options.preproc_prefix, '.nii']);
for it = 1:length(listing_s)
    tmp_file =listing_s(it).name;
    if ~isempty(tmp_file)
        count = count + 1;
        raw_filenames{count} = [nifti_dir, session_dir, '/', tmp_file]; %#ok<AGROW>
    end
end

% read spm volume using file names
subj = load_spm_pattern(subj, 'epi', [roi_type '_' roi_name], raw_filenames);



%% save subj structure

if ~isempty(options.save_dir)
    save(options.save_dir, 'subj')
end



end


