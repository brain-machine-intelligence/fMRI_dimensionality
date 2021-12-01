function [boldpat_vet, N_trial] = arbMBMF_boldpat(Exp, ROI, ID, varargin)
% return EPI volume as masked time series

% boldpat_vet: roi_size x length(event) x N_trial
% multi-dimensional array with Voxel, Event, Trial dimension

% =========================================================================
% Name-value pairs
% 'TrialCondition' - 
% 'Event' - within-trial event index (1:f1, 2:S1, 3:A1, 4:f2, 5:S2, 6:A2, 7:f3, 8:R)
%           possible range: [-7:0, 1:8, 9:16] (prev, curr, next trial)
% 'roi_type' - 'sphere10', 'sphere5', 'AAL3'

% =========================================================================

% root = '/home/ydsung'; is_cluster = 1;
root = '//143.248.30.94/bmlsamba/ydsung'; is_cluster = 0;

% default input parameters
options = struct('SignalType', 'zscore', ...
                'MaskNum', [], ...
                'Event', 1:8, ...
                'TR', 2.78, ...
                'Delay', 6.1, ...
                'TrialCondition', 'none', ...
                'strict', [], ...
                'interpolation', [], ...
                'roi_type', 'sphere10', ...
                'fextension', '.nii', ...
                'prefix', 'wra*');
% read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(arbMBMF_boldpat) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(arbMBMF_boldpat) %s is not a recognized parameter name', pair{1})
    end
end

%% Experiment type

switch Exp
    
    case 'Lee2014'
        % Lee2014 subject
        % sbj ids
        IDs = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};  % 20:two small sessions(2)

        % SBJstructure
%         sbj_struct = load([root '/A_Research/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
        sbj_struct = ... 
            load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\SBJ_structure_each_exp_BESTcollectionJan29');
        SBJ = sbj_struct.SBJ{ID};
        
        % Session : 1~2, 1~3, 1~4 or 1~5
        N_sessions = 5*ones(1,24); N_sessions(11)=4; N_sessions(16)=3; N_sessions(20)=2; N_sessions(21)=4;
        N_session = N_sessions(IDs{ID});
        % Session marks
        sess_mark_struct = ... 
            load([root '/A_Research/Dim_control_PFC_metaRL/session_marks/sess_mark_Lee2014.mat']);     % load 'sess_mark_Lee2014'
        sess_mark_Lee2014 = sess_mark_struct.sess_mark_Lee2014;
        sess_marks = sess_mark_Lee2014(IDs{ID},:);
        
        % neuroimaging file (functional EPI volume) directory
        if strcmp(root, '/home/ydsung')
            fMRI_dir_format = ...
                [root '/2014fmri/fmri_arbitration/od-arbitration-%03d/func/run_%04d/'];
        else
            fMRI_dir_format = 'D:/fmri_arbitration/od-arbitration-%03d/func/run_%04d/';
        end
        
    case 'Heo2018'
        
        IDs = {1, 3, 11, 12, 18, 19, 20, 24, 25, 27, 32, 33, 34, 35, 37, 38, 40, 41, 42, 44, 46, 48, 50, 51, 52, 53, 54, 55};
        % SBJstructure
%         sbj_struct = load([root '/A_Research/Dim_control_PFC_metaRL/SBJ_structure']);
%         sbj_struct = sbj_struct.SBJ2;
        sbj_struct = load([root '/A_Research/Dim_control_PFC_metaRL/SBJ_structure_1027_for_revision']);
        sbj_struct = sbj_struct.SBJ3;
        SBJ = sbj_struct{IDs{ID}};
        
        % Session : 1~4
        N_sessions = 4*ones(1,55);
        N_session = N_sessions(IDs{ID});
        % Session marks
        sess_mark_struct = ... 
            load([root '/A_Research/Dim_control_PFC_metaRL/session_marks/mark_depression.mat']);     % load 'sess_mark_Lee2014'
        mark_depression = sess_mark_struct.mark_depression;
        sess_marks = mark_depression(IDs{ID},:);
        
        % neuroimaging file (functional EPI volume) directory
        if strcmp(root, '/home/ydsung')
            fMRI_dir_format = ...
                '/home/syh1743/Depression_ana/Sbj%d/od-arbitration-001/functional/session%d/';
        else
            fMRI_dir_format = ...
                'D:/fmri_depression_ana/Sbj%d/session%d/';
        end
        
                
    case 'Kim2019'
        
        IDs = 1:24;  IDs([9, 12, 16]) = [];
        IDs = num2cell(IDs);
                
        if strcmp(root, '/home/ydsung')
%             root2 = '/home/bmlshare';
            root2 = '/home/ydsung';   % 2021-05-05
            sbj_struct = load([root2 '/complexity/modelRLsource/result_simul/SBJ_structure.mat']);
        else
%             root2 = '//143.248.30.94/bmlsamba/bmlshare';
%             sbj_struct = load([root2 '/complexity/modelRLsource/result_simul/SBJ_structure.mat']);
%             root2 = '//143.248.30.94/bmlsamba/bmlshare';
            root2 = '//143.248.30.94/bmlsamba/ydsung';    % 2021-05-05
            sbj_struct = ...
                load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_Kim2019.mat');
        end
                
        % SBJ structure
        SBJ = sbj_struct.SBJ{ID};
        
        % neuroimaging file (functional EPI volume) directory
        fMRI_dir_format = ...
            [root2 '/complexity/od-cog-%03d/func/run_%04d/'];
        
        % Session
        N_sessions = [5,5,4,6,5,5,6,6,-1,5,5,-1,5,5,5,-1,6,4,5,6,5,5,5,5];
        N_session = length(dir(sprintf([root2 '/complexity/od-cog-%03d/func/'], ...
            IDs{ID}))) - 2;
        % sanity check
        if N_session ~= N_sessions(IDs{ID}); error('(arbMBMF_boldpat) N_session error'); end
        
        % Session marks
        mark_path = [root sprintf('/A_Research/Dim_control_PFC_metaRL/session_marks/mark_Kim2019_%d.mat', ...
            IDs{ID})];
        try            
            sess_mark_struct = load(mark_path);
            sess_marks = sess_mark_struct.mark_Kim2019;
        catch
            disp('(arbMBMF_boldpat) mark (# of volumes of each run) counting')
            mark_Kim2019 = zeros(1, N_session+1);   % initialization
            for ss = 1:N_session
                n_files = dir([sprintf(fMRI_dir_format, IDs{ID}, ss), options.prefix, '.nii']);
                mark_Kim2019(ss+1) = length(n_files);
            end
            mark_Kim2019 = cumsum(mark_Kim2019);
            save(mark_path, 'mark_Kim2019');
            sess_marks = mark_Kim2019;
        end
        
    case 'Shin2019'
        
        IDs = [5:8,11,13,15:26];
        IDs = num2cell(IDs);
%         od_list=[5:8,11,13,15:26]
        
        if strcmp(root, '/home/ydsung')
            root2 = '/home/SJH_indi';
            sbj_struct = load([root '/A_Research/Dim_control_PFC_metaRL/SBJ_struture_Shin2019.mat']);
        else
            root2 = '\\143.248.30.94\bmlsamba\SJH_indi';
            sbj_struct = ...
                load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_struture_Shin2019.mat');
        end
        
        % SBJ structure
        SBJ = sbj_struct.SBJ{ID};
        
        % neuroimaging file (functional EPI volume) directory
        fMRI_dir_format = ...
            [root2 '/Human_Guidance/fMRI/images/od-arbitration-%03d/func/run_%04d/'];
        
        % Session
%         N_sessions = [5,5,4,6,5,5,6,6,-1,5,5,-1,5,5,5,-1,6,4,5,6,5,5,5,5];
        N_session = length(dir(sprintf([root2 '/Human_Guidance/fMRI/images/od-arbitration-%03d/func/'], ...
            IDs{ID}))) - 2;
        % sanity check
%         if N_session ~= N_sessions(IDs{ID}); error('(arbMBMF_boldpat) N_session error'); end
        
        % Session marks
        mark_path = [root sprintf('/A_Research/Dim_control_PFC_metaRL/session_marks/mark_Shin2019_%d.mat', ...
            IDs{ID})];
        try            
            sess_mark_struct = load(mark_path);
            sess_marks = sess_mark_struct.mark_Shin2019;
        catch
            disp('(arbMBMF_boldpat) mark (# of volumes of each run) counting')
            mark_Shin2019 = zeros(1, N_session+1);   % initialization
            for ss = 1:N_session
                n_files = dir([sprintf(fMRI_dir_format, IDs{ID}, ss), options.prefix, '.nii']);
                mark_Shin2019(ss+1) = length(n_files);
            end
            mark_Shin2019 = cumsum(mark_Shin2019);
            save(mark_path, 'mark_Shin2019');
            sess_marks = mark_Shin2019;
        end
        
end

% ROI
% ROIs: {'FPC', 'lilPFC', 'rilPFC', 'rACC', 'pPut', 'omPFC', 'vmPFC', 'lV1', 'rV1', 'lIT', 'rIT', 'lCaudate', 'rCaudate', 'lHPC', 'rHPC'};
% ROI = 'lilPFC'; % disp(['sbj' num2str(ID) ' ' ROI]);

%% BOLDpattern

switch options.SignalType
    case 'zscore'
        boldpat_dir = [root '/A_Research/Dim_control_PFC_metaRL/subj_masked_EPI/' ... 
            Exp '/' options.roi_type '/' erase(options.prefix, '*') '/' ROI];
        
        if strcmp(options.roi_type, 'AAL3')
            fname_format = ['subj_EPI_' Exp '_%s_%s'];
        else
            fname_format = 'subj_fMRI_sbj14md14mk14_%s_%s';
        end
        
        % load BOLD pattern
        %     boldpat_name = [boldpat_dir '/subj_fMRI_sbj14md14mk14_' ROI '_' num2str(IDs{ID}) '.mat'];
        % boldpat_name = [boldpat_dir '/' sprintf(fname_format, ROI, num2str(IDs{ID})+id_jump) '.mat'];
        boldpat_name = [boldpat_dir '/' sprintf(fname_format, ROI, num2str(IDs{ID})) '.mat'];
        try
            subj_struct = load(boldpat_name);     % try to load 'subj' struct
            subj = subj_struct.subj;
            boldpat = get_mat(subj,'pattern','epi_z');
        catch
            % Activated only for the first run & use saved data in later runs
            disp(['---------- EPI pattern loading for ' ROI ' ----------'])
%             % subj_init_with_fMRI(exp_name, region, sbj_id, N_session, maskDir, maskType, fMRIdir_format, fMRIprefix, saveDir)
%             [subj, ~] = subj_init_with_fMRI(Exp, ROI, IDs{ID}, N_session, ...
%                 [root '/A_Research/princeton-mvpa-toolbox-master/ROI_mask'], ... % maskDir
%                 options.roi_type, ...           % maskType
%                 fMRI_dir_format, ...            % funcitonal EPI volume directory format
%                 options.prefix, ...             % preproc prefix
%                 boldpat_name);                  % fname to save 'subj' struct
%             boldpat = get_mat(subj,'pattern','epi_z');

            % 2021-03-20
            boldpat = [];
            for is = 1:N_session
                
                % For the new mask
                % subj = roi_nifti_load_2(exp, sbj_id, session, roi_type, roi_name, varargin)
                subj = roi_nifti_load_2(Exp, ID, is, options.roi_type, ROI, ...
                    'is_cluster', is_cluster, ...
                    'preproc_prefix', options.prefix, ...
                    'mask_num', options.MaskNum, ...
                    'save_dir', []);                

                sel_name = sprintf('sel%d', is);
                which_session = ones(1, size(get_mat(subj,'pattern','epi'), 2));
                subj = initset_object(subj, 'selector', ...
                    sel_name, which_session);
                
                % Detrending
                subj = detrend_pattern(subj, 'epi', sel_name);
                
                % Z scoring
                subj = apply_to_runs(subj, 'epi_dt', sel_name, 'apply_zscore', ... 
                    'new_patname', sprintf('epi_z_%d', is));
                
                ts_z = get_mat(subj, 'pattern', sprintf('epi_z_%d', is));
                boldpat = [boldpat, ts_z];

            end % session
            subj = init_object(subj,'pattern', 'epi_z');
            subj = set_mat(subj, 'pattern', 'epi_z', boldpat);
            save(boldpat_name, 'subj');
            
        end
        % boldpat = get_mat(subj,'pattern','epi');
        % boldpat = get_mat(subj,'pattern','epi_dt');
        % boldpat = get_mat(subj,'pattern','epi_z');
        
    case 'percent'
        boldpat = [];
        for is = 1:N_session
            subj_save_path = ['Z:\JR\dimensionality/subj/' options.roi_type '/' ROI '/subj_' num2str(ID) '_' num2str(is)];
            try
                subj_struct = load(subj_save_path);
                subj = subj_struct.subj;
            catch
                % For the new mask
                % subj = roi_nifti_load_2(exp, sbj_id, session, roi_type, roi_name, varargin)
                subj = roi_nifti_load_2(Exp, ID, is, options.roi_type, ROI, ...
                    'is_cluster', is_cluster, ...
                    'mask_num', options.MaskNum, ...
                    'save_dir', subj_save_path);
            end
            
            % get percent signal
            ts = get_mat(subj,'pattern','epi');
            ts_p = (ts - mean(ts,2)) ./ (mean(ts,2) * ones(1,size(ts,2)));
            
            boldpat = [boldpat, ts_p];
            
        end % session
        
        % voxel screening
        var_ts_p = var(boldpat, [], 2);
        ts_p_thr = prctile(var_ts_p, 95);
        boldpat = boldpat(var_ts_p < ts_p_thr, :) * 100;
        
end


%% Finding fMRI vol index of cue1/cue2/reward/reward_end in each trial (based on HIST_event_info)
TR = options.TR;
delay = options.Delay;

MRI_ids = [];
blk_con = [];
for sess = 1:N_session
%     mri_ids = sess_marks(sess) + SARtime_to_MRIidx(SBJ.HIST_event_info{sess}, TR, delay);
    mri_ids = sess_marks(sess) + ... 
        SARtime_to_MRIidx(SBJ.HIST_event_info{sess}, TR, delay, ... 
        'Exp', Exp, ...
        'strict', options.strict, ...
        'forInterpol', options.interpolation);
    mri_ids(mri_ids>sess_marks(sess+1)) = sess_marks(sess+1);
    MRI_ids = [MRI_ids, mri_ids]; % 8 x N_trials
    blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
end
% 'MRI_ids' expansion (add previous & next trial)
MRI_prev = MRI_ids; MRI_next = MRI_ids;
MRI_prev(:,2:end) = MRI_ids(:,1:(end-1)); MRI_prev(:,1) = nan;
MRI_next(:,1:(end-1)) = MRI_ids(:,2:end); MRI_next(:,end) = nan;
MRI_ids = [MRI_prev; MRI_ids; MRI_next]; % 24 x N_trials

switch options.TrialCondition
    case 'none'
        N_trial = size(MRI_ids,2);
        tri = 1:N_trial;
    case 1 % 1:specific goal, 0:flexible goal
        tri = ismember(blk_con, [1,2]);
    case 2 % 1:flexible goal, 0:specific goal
        tri = ismember(blk_con, [3,4]);
    case 3 % 1:low uncertainty, 0:high uncertainty
        tri = ismember(blk_con, [1,4]);
    case 4 % 1:high uncertainty, 0:low uncertainty
        tri = ismember(blk_con, [2,3]);
end

% possible input events: -7~16, corresponding idx: 1~24
event = options.Event;

if iscell(event)
    % example: event = {[2 3 4], [5 6 7],[8 9]};  
    boldpat_vet = nan * ones(size(boldpat,1), length(event), ...
        size(MRI_ids(:,tri), 2));
    for ei = 1:length(event)
        full_idx = reshape(MRI_ids(event{ei}+8, tri), 1, []);
        null_idx = isnan(full_idx);
        new_idx = full_idx; new_idx(null_idx) = 1;
        boldpat_temp = boldpat(:, new_idx);
        boldpat_temp(:, null_idx) = nan;
        boldpat_temp = reshape(boldpat_temp, ...
        size(boldpat,1), length(event{ei}), []);
%         boldpat_vet(:, ei, :) = mean(boldpat_temp, 2);
        boldpat_vet(:, ei, :) = nanmean(boldpat_temp, 2);
    end
else
    % example: event = 1:8; event = [-1 3 15];  
    % multi-dimensional BOLD pattern (N_voxel x N_event x N_trial)
    full_idx = reshape(MRI_ids(event+8, tri), 1, []);
    null_idx = isnan(full_idx);
    new_idx = full_idx; new_idx(null_idx) = 1; % arbitrary positive integer
    boldpat_vet = boldpat(:, new_idx);
    boldpat_vet(:, null_idx) = nan;
    boldpat_vet = squeeze(reshape(boldpat_vet, ...
        size(boldpat,1), length(event), []));
    N_trial = size(boldpat_vet, 3);
    % disp(size(boldpat_event)); % check
    % boldpat_event = permute(boldpat_event, [1 3 2]); % 200602
end


%% Supplementary info

% Subject =================================================================
% === Lee2014 ===
% SBJ_NAME = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark', ...
%     'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
% SBJ_ID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};
% === Heo2018 ===
% SBJ_NAME = {'subject001_ksy',[],'subject03_keb','subject004_ksj','subject5_lwc','subject6_sjh','subject7_kch','subject8_kjs', ...
%     'subject9_lks','subject10_ssh','subject11_hjy','subject12_khs','subject13_pjb','subject14_syh','subject15_kik','subject16_jja', ...
%     'subject16_lsh','subj18_kjh','subj19_kny','subj20_jjh','subject21_lsl','subject22_jsh','subject23_yhj','subj24_kjy','subj25_kjs', ...
%     'subject26_cjb','subj27_kkm','subject28_ljh','subject29_cyj','subject30','subj31_ssw','subj32_khj','subj33_ohy','subj34_yhw','subj35_ajs', ...
%     [],'subj37_shi','subj38_lsh',[],'subj40_kkr','subj41_jhj', 'sbj42_ljk',[],'sbj44_ksy', 'sbj45_jej','sbj46_ses','sbj47_hej','sbj48_kdy', ...
%     'subject49_ljy','sbj50_kdh','sbj51_akr', 'sbj52_jkw','sbj53_nyk','sbj54_jyh', 'sbj55_ker'};
% SBJ_ID = {1, 3, 11, 12, 18, 19, 20, 24, 25, 27, 32, 33, 34, 35, 37, 38, 40, 41, 42, 44, 46, 48, 50, 51, 52, 53, 54, 55};
% === Kim2019 ===
% 

% ROI =====================================================================
% Insula & lPFC (SPE)
% v.Striatum & d.Striatum-putamen (RPE)
% ilPFC (RelMB,RelMF,max)
% p.Put (QMF)
% omPFC (QMB)
% vmPFC (Qarb)
% === 5mm sphere ===
% ROI = {'FPC', 'lilPFC', 'lV1', 'omPFC', 'pPut', 'rACC', 'rilPFC', 'rV1', 'vmPFC'};
% === 10mm sphere ===
% ROI = {'FPC', 'lilPFC', 'rilPFC', 'rACC', 'pPut', 'omPFC', 'vmPFC', 'lV1', 'rV1', 'lIT', 'rIT', 'lCaudate', 'rCaudate', 'lHPC', 'rHPC'};

% Trial condition (within-subject effect) =================================
% === Block condition ===
% Goal condition (specific vs flexible), Uncertainty (low vs high),
% Complexity (low vs high)
% === Behavioral performance ===
% Hit (rewarded vs no rewarded), Optimality (optimal vs non-optimal)
% === Learning stage? (long-term temporal change) ===
% early vs late, each session
% === Value update amount? ( SIGMA_s,a(Q(s,a)) ) ===
% large update vs smalle update


