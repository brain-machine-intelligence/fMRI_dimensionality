
function [] = arbMBMF_shattering(id)

% addpath(genpath('/home/ydsung/A_Research'));

% Setting =================================================================

% ===== experiment =====
exp = 'Lee2014';  exp_sfx = [];       prefix = 'wra*';
% exp = 'Heo2018';
% exp = 'Kim2019';    exp_sfx = '_Kim_wa';    prefix = 'wa*';
% exp = 'Shin2019';


% ===== ROI mask for multi-voxel pattern (MVP) =====
% ROI = {'ilPFC', 'V1'};
% ROI = {'lilPFC'};
% ROI = {'lilPFC', 'FPC', 'vmPFC', 'lV1'};
% ROI = {'rilPFC', 'rV1', 'vmPFC'};
% ROI = {'rilPFC', 'FPC'};
% ROI = {'lilPFC', 'rilPFC', 'FPC', 'vmPFC', 'lV1', 'rV1'};
% ROI = {'lilPFC', 'rilPFC', 'FPC', 'vmPFC', 'lV1', 'rV1', 'lIT', 'rIT'};

% Anatomical ROI (AAL3)
% ROI = {'F3TL'};  % for test 201020
% ROI = {'V1L'};  % for test 201020
% ROI = {'F3TL', 'F3TR', 'G_FIT'};
% ROI = {'V1L', 'V1R'};
% ROI = {'lilPFC', 'rilPFC', 'F3TL', 'F3TR', ... 
%     'V1L', 'V1R', 'VentStrL', 'VentStrR', 'HIPPOL', 'HIPPOR'};
% ROI = {'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};
% ROI = {'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'};

% +DLPFC & R/L integrated OFC
ROI = {'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', 'DLPFCL', 'DLPFCR', 'OFC'};
% ROI = {'DLPFCL', 'DLPFCR', 'OFC'};

% 201029
% ROI = {'V1L', 'V1R'}; % sphere10
% ROI = {'lilPFC', 'rilPFC', 'F3TL', 'F3TR', 'V1L', 'V1R', 'V1L', 'V1R'};
% ROI = {'F3TL'};
% ROI = {'F3TR'};

% 201102
% ROI = {'F3TR'};
% ROI = {'F3TL', 'F3TR'};
% ROI = {'VentStrL', 'VentStrR'};

% 201110
% ROI = {'PALLL', 'PALLR', 'PutL', 'PutR', 'CdL', 'CdR'};

% % 201111
% rng(2020)
% N_us = 17; % voxel undersample (without replacing)
% % N_us = 15; % voxel undersample (without replacing)
% % N_us = 6; % voxel undersample (without replacing)
% N_vox = 2529;
% % N_vox = 2258;
% % N_vox = 932;
% N_min_vox = 143;
% rand_vox_ids = randperm(N_vox);
% sub_vox = cell(1, N_us);   % subsets of voxel index
% ROI = cell(1, N_us);
% roi_name = 'F3TL';
% % roi_name = 'V1L';
% % roi_name = 'HIPPOL';
% for ius = 1:N_us
%     ROI{ius} = [roi_name 'us' num2str(ius)];
%     sub_vox{ius} = rand_vox_ids((ius-1)*N_min_vox+1:ius*N_min_vox);
% end



% ====================== labelling by task variable =======================
% CClabel_name = 'UncCond';
% CClabel_name = 'CmplxCond';
% CClabel_name = 'Goal(6,7,8)';
CClabel_name = [];

% LABEL = {'R(0,20,40)_G uH', 'R(0,20,40)_H uH'};     % 2021-04-17 fix the uncCond: uH
% LABEL = {'Session'};    % 2021-04-15 fix the uncCond: uH
% LABEL = {'S(6,7,9)_G uH', 'S(6,7,9)_H uH'};     % 2021-04-15 fix the uncCond: uH

% LABEL = {'R(0,20,40)_G', 'R(0,20,40)_H'};     % 2021-04-09 

% LABEL = {'blkCond early', 'blkCond late'};     % 2021-04-08 
% LABEL = {'R(0,20,40) binPMB 0', 'R(0,20,40) binPMB 1'};     % 2021-04-08

% LABEL = {'S3(6,7,9) binPMB perC 0', 'S3(6,7,9) binPMB perC 1', ...     % 2021-04-06
%          'R(0,20,40) binPMB perC 0', 'R(0,20,40) binPMB perC 1', ...   % 2021-04-06
%          'R binPMB perC 0', 'R binPMB perC 1'};                        % 2021-04-06

% LABEL = {'Coin(1,2,3) binPMB 0', 'Coin(1,2,3) binPMB 1'};   % 210401
% LABEL = {'S3(6,7,9) binPMB 0', 'S3(6,7,9) binPMB 1'};   % 210401
% LABEL = {'S3(6,7,9) binChoOpt2 perC 1', 'S3(6,7,9) binChoOpt2 perC 0'};    % 210401

% LABEL = {'Goal(1,2,3) x Unc'};     % 210330 for Kim2019 
% LABEL = {'Goal(1,2,3) x Cmplx'};     % 210330 for Kim2019 

% LABEL = {'Unc ChoOpt 0', 'Unc ChoOpt 1'};     % 210318
% LABEL = {'Unc in H'};     % 210318

% LABEL = {'Goal(6,7,8) x Unc'};
LABEL = {'Goal(6,7,8) x Unc LOROV'};

% LABEL = {'S3G(6,7,8)', 'RG(10,20,40)'};     % 210312
% LABEL = {'S3 new'};     % 210312

% LABEL = {'Goal(6,7,8) uL', 'Goal(6,7,8) uH'};     % 210311
% LABEL = {'Goal(6,7,8)'};     % 210311
% LABEL = {'R'};     % 210311
% LABEL = {'Goal uL', 'Goal uH'};   % 210311

% ===== task-relevant context mixing =====
% LABEL = {'Goal x Unc'};

% ===== temporal realignment of context embedding =====
% LABEL = {'CmplxCond P1', 'CmplxCond P2', 'CmplxCond P3', ...
%         'UncCond P1', 'UncCond P2', 'UncCond P3', ...
%         'Interaction P1', 'Interaction P2', 'Interaction P3'};    % 2021-05-14 Kim2019
% LABEL = {'GoalCond P1', 'GoalCond P2', 'GoalCond P3', ...
%         'UncCond P1', 'UncCond P2', 'UncCond P3', ...
%         'Interaction P1', 'Interaction P2', 'Interaction P3'};    % 2021-05-12 Lee2014

% ===== basic taskVar =====
% LABEL = {'S3 run2345'};     % 2021-06-03 label runs test
% LABEL = {'GoalCond', 'UncCond', 'Interaction'};    % 2021-05-09 Lee2014
% LABEL = {'UncCond'};
% LABEL = {'GoalCond LOROV', 'UncCond LOROV', 'Interaction LOROV'};
% LABEL = {'CmplxCond', 'UncCond', 'Interaction'};  % 2021-05-09 Kim2019
% LABEL = {'blkCond'};
% LABEL = {'S2'};
% LABEL = {'S3'};
% LABEL = {'blkCond', 'S2', 'S3'};
% LABEL = {'blkCond', 'S2', 'Coin'};  % 2021-03-25
% LABEL = {'Goal'};
% LABEL = {'Goal(6,7,8) LOROV', 'Goal(6,7,8) uL LOROV', 'Goal(6,7,8) uH LOROV'};
% LABEL = {'Goal(6,7,8) LOROV'};
% LABEL = {'Goal(6,7,8) uL LOROV', 'Goal(6,7,8) uH LOROV'};
% LABEL = {'Goal cL LOROV', 'Goal cH LOROV', 'Goal LOROV'};
% LABEL = {'Goal LOROV'};
% LABEL = {'Goal uL LOROV', 'Goal uH LOROV'};
% LABEL = {'Goal uL LOROV', 'Goal uH LOROV', 'Goal cL LOROV', 'Goal cH LOROV'};
% LABEL = {'S2', 'blkCond'};
% LABEL = {'blkCond', 'S3'};
% LABEL = {'S3', 'blkCond'};
% LABEL = {'S3 shuff', 'blkCond shuff'};

% ===== latent RL signal =====
% LABEL = {'RPE', 'SPE'};     % 2021-05-22
% LABEL = {'binRelMB', 'binRelMF', 'binRelMAX'};     % 2021-06-02

% ==== BHV related embedding ====
% 2021-11-15 sbj1~6 run finish
% LABEL = {'Goal(6,7,8) binPMB_G0 LOROV', 'Goal(6,7,8) binPMB_G1 LOROV'}; 
% LABEL = {'Goal(6.7.8) binChoOpt_G0 LOROV', 'Goal(6.7.8) binChoOpt_G1 LOROV'}; 
% LABEL = {'Goal(6,8) binChoOpt_G0 LOROV', 'Goal(6,8) binChoOpt_G1 LOROV'}; 

% ===== validation =====
% LABEL = {'S3() after S3 7'};    % 201212, hyp: abstract context info vs context-dependent S3 bias effect
% LABEL = {'blkCond in S2 4'};    % 201212, hyp: abstract context info vs context-dependent S2 bias effect
% LABEL = {'blkCond in A2 right'};    % 201212, hyp: abstract context info vs context-dependent action bias effect
% LABEL = {'blkCond only goal 6,-1'};    % 201212, hyp: abstract context info vs context-dependent goal bias effect

% LABEL = {'blkCond in S3 7'};    % 201212, hyp: abstract context info vs context-dependent S3 bias effect
% LABEL = {'blkCond in S2 4'};    % 201212, hyp: abstract context info vs context-dependent S2 bias effect
% LABEL = {'blkCond in A2 right'};    % 201212, hyp: abstract context info vs context-dependent action bias effect
% LABEL = {'blkCond only goal 6,-1'};    % 201212, hyp: abstract context info vs context-dependent goal bias effect

% 201203 event 1 2 4 5 7 9 10
% LABEL = {'S3(6,7,9) binQarb2 0', 'S3(6,7,9) binQarb2 1', ... 
%     'blkCond_binrelMAX_0', 'blkCond_binrelMAX_1'};

% 201125
% LABEL = {'S3(6,7,9) binQarb2 0', 'S3(6,7,9) binQarb2 1', ... 
%     'blkCond_binChoOpt0', 'blkCond_binChoOpt1'};

% 201124
% LABEL = {'blkCond_binrelMAX_0', 'blkCond_binrelMAX_1'};

% 201123
% LABEL = {'S3 binChoOpt2 perC 0', 'S3 binChoOpt2 perC 1'};
% LABEL = {'blkCond binChoOpt2 perC 0', 'blkCond binChoOpt2 perC 1'};   % 201208

% 201110 - BG ROIs
% LABEL = {'S3', 'blkCond'};
% LABEL = {'S3'};
% LABEL = {'blkCond'};

% 201104
% LABEL = {'blkCond ChoCon 0', 'blkCond ChoCon 1'};

% 201103
% LABEL = {'binrelMAX3xrelComp3'};
% LABEL = {'blkCond_binrelMAX_0', 'blkCond_binrelMAX_1'};

% 201102
% LABEL = {'S2'};
% LABEL = {'S3'};     % for event -1 0 9 10

% 201029
% LABEL = {'blkCond shuff'};
% LABEL = {'S3'};
% LABEL = {'blkCond'};

% % for test 201020
% LABEL = {'S3'};

% 201020 ~ 201027? slave1
% LABEL = {'S3', 'S3 shuff'};


% LABEL = {'currBlk blkchng shuff', 'currBlk blkchng+1T shuff', 'currBlk blkchng+2T shuff', ...
%     'currBlk blkchng+3T shuff', 'currBlk blkchng+4T shuff'};
% LABEL = {'prevBlk blkchng shuff', 'prevBlk blkchng+1T shuff', 'prevBlk blkchng+2T shuff', ...
%     'prevBlk blkchng+3T shuff', 'prevBlk blkchng+4T shuff'};

% LABEL = {'currBlk blkchng', 'currBlk blkchng+1T', 'currBlk blkchng+2T', 'currBlk blkchng+3T', ...
%     'prevBlk blkchng', 'prevBlk blkchng+1T', 'prevBlk blkchng+2T', 'prevBlk blkchng+3T'};
% LABEL = {'currBlk blkchng', 'currBlk blkchng+1T', 'currBlk blkchng+2T', ...
%     'currBlk blkchng+3T', 'currBlk blkchng+4T'};
% % LABEL = {'prevBlk blkchng', 'prevBlk blkchng+1T', 'prevBlk blkchng+2T', ...
% %     'prevBlk blkchng+3T', 'prevBlk blkchng+4T'};


% LABEL = {'S3(6,7,9) binQarb2 0', 'S3(6,7,9) binQarb2 1'};

% LABEL = {'blkCond prevT', 'blkCond nextT'};

% LABEL = {'blkCond early', 'blkCond late'};
% LABEL = {'binrelMAX3xrelComp3 early', 'binrelMAX3xrelComp3 late'};

% LABEL = {'S3(6,7,9)_G new', 'S3(6,7,9)_H new'};
% LABEL = {'S3(6,7,9)_G', 'S3(6,7,9)_H'};
% LABEL = {'S3(6,7,9)_uH', 'S3(6,7,8)_uL'}; 

% LABEL = {'S3 early shuff', 'S3 late shuff'};
% LABEL = {'S3 early', 'S3 late'};
% LABEL = {'S3(6,7,9)', 'S3(6,7,9) shuff'};
% LABEL = {'S3(6,7,9)_G', 'S3(6,7,9)_G shuff'};

% LABEL = {'S3G(6,7,8)_uL shuff', 'S3G(6,7,8)_uH shuff'};
% LABEL = {'S3G(6,7,8)_uL', 'S3G(6,7,8)_uH'}; % S3G per UncCond
% LABEL = {'S3(6,7,9)_H shuff', 'S3(6,7,9)_uH shuff', 'S3(6,7,8)_uL shuff'};
% LABEL = {'S3(6,7,9)_H', 'S3(6,7,9)_uH', 'S3(6,7,8)_uL'}; % S3 per GoalCond/UncCond
% LABEL = {'S3'};

% LABEL = {'blkCond_binrelMAX_0', 'blkCond_binrelMAX_1'};
% LABEL = {'S3(6,7,9)_binChoOpt2_0', 'S3(6,7,9)_binChoOpt2_1', ...
%     'S3(6,7,9)_binChoOpt_0', 'S3(6,7,9)_binChoOpt_1'};
% LABEL = {'blkCond_binChoOpt0', 'blkCond_binChoOpt1'};
% LABEL = {'blkCond', 'S2G', 'S3G', 'A1', 'A2', 'Hit', 'ChoOpt1', 'ChoOpt2'};
% LABEL = {'S3', 'blkCond'};
% LABEL = {'S3 shuff', 'S3G shuff', 'blkCond shuff'};
% LABEL = {'S3G', 'blkCond', 'binrelMB3xbinrelMF3'};
% LABEL = {'S3G', 'blkCond', 'binrelMAX3xrelComp3', 'binrelMB3xbinrelMF3'};
% LABEL = {'A1', 'A2'};
% LABEL = {'S3_blk1', 'S3_blk2'};
% LABEL = {'Hit'};
% LABEL = {'blkCond'};
% LABEL = {'GoalCond'};
% LABEL = {'Session'};
% LABEL = {'blkCond_sess1', 'blkCond_sess2', 'blkCond_sess3', 'blkCond_sess4', 'blkCond_sess5'};
% LABEL = {'relMB2', 'relMB3', 'relMF2', 'relMF3', ...
%     'relMAX2', 'relMAX3', 'GoalCond', 'UncCond'};
% LABEL = {'binrelMAX3xrelComp3'};
% LABEL = {'binrelMAX3xrelComp3_sess1', 'binrelMAX3xrelComp3_sess2', ... 
%     'binrelMAX3xrelComp3_sess3', 'binrelMAX3xrelComp3_sess4', ... 
%     'binrelMAX3xrelComp3_sess5'};


% N_allclass = [3 4];
% N_allclass = 4;
% N_allclass = [4 2 4];
% N_allclass = [3 3 4 4];
% N_allclass = [3 3 3 3 4 4];
% N_allclass = 2 * ones(1, length(LABEL));
% N_allclass = 3 * ones(1, length(LABEL));
% N_allclass = 5 * ones(1, length(LABEL));
% N_allclass = 4 * ones(1, length(LABEL));
N_allclass = 6 * ones(1, length(LABEL));
% N_allclass = 8 * ones(1, length(LABEL));

% to categorize continuous variables into several classes
categorize = cell(1, length(LABEL));  % empty cell for no categorization
% categorize = num2cell(2 * ones(1, length(LABEL)));
% categorize = num2cell(3 * ones(1, length(LABEL)));
% categorize = num2cell(N_allclass);

% lab_shuffle = 1;
lab_shuffle = 0;



% ========================= MVP setting =========================
% prefix = 'wra*';  % preprocessing prefix
% prefix = 'swar*';
% prefix = 'swa*';
% prefix = 'wa*';

% SignalType = 'percent';
SignalType = 'zscore'; % default setting

% ROI_types = {};
% ROI_types = {'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3', 'AICHAmc'};
% ROI_types = {'sphere10', 'sphere10', ..., 
%     'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
% ROI_types = {'sphere10', 'sphere10', 'AAL3', 'AAL3', ..., 
%     'sphere10', 'sphere10', 'AAL3', 'AAL3'};
% ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3'};
ROI_types = {'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL3', 'AAL', 'AAL', 'AAL'}; % +DLPFC & R/L integrated OFC
% roi_type = 'sphere5';
% roi_type = 'sphere10';
% roi_type = 'AAL3';
% roi_type = 'AAL';
roi_type = [];

if ~isempty(roi_type)
    switch roi_type
        case 'AAL3'
            MaskNum_map = containers.Map({'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR'}, ...
                [9       10     47     48      41         42        157        158]);
        case 'AAL'
            MaskNum_map = containers.Map({'DLPFCL', 'DLPFCR', 'OFCL', 'OFCR', 'OFC'}, ...
                {7  8   [5 9 15 27]  [6 10 16 28]  [5 9 15 27 6 10 16 28]});
    end
else
    MaskNum_map = containers.Map({'F3TL', 'F3TR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', ...
        'DLPFCL', 'DLPFCR', 'OFC'}, ...
        {9       10     47     48      41         42        157        158  ...
        7  8   [5 9 15 27 6 10 16 28]});
end

% EVENT = 5;
% EVENT = 8;
% EVENT = [-4 13];        % for test 201020
% EVENT = {[1 2], [4 5], [7 8]};
% EVENT = {[1 2 4 5 7 8]};
% EVENT = [-1 0];                 % make-up run for F3TR S3
% EVENT = [9 10];               % make-up run for F3TL & F3TR S3
% EVENT = [1 2 4 5 7 8];        % mainly for 'strict' MVP
% EVENT = [1 2 4 5 7 9 10];        
EVENT = [1 2 4 5 7 8 9 10];        % maing setting
% EVENT = [-1 0 1 2 4 5 7 8 9 10];
% EVENT = [7 8 9 10];     % make-up run for F3TL blkCond (slave1 too slow)
% EVENT = -3:8;
% EVENT = -3:12;
% EVENT = 1:8;
% EVENT = {[2 3], [5 6], 8};

strict = 0;
% strict = 1;
interpolation = 0;



% ========================= fitclinear setting =========================
% Repeat = 1;    % # of undersample
Repeat = 20;    % # of undersample
% Repeat = 50;    % # of undersample

% rRepeat = 10;
rRepeat = 0;

% Coding = 'onevsone';
% Learners = 'svm';
% KFold = 10;
% KFold = 5;
% KFold = 3;
KFold = [];

% CVtype = 'LOOV';
CVtype = 'LOROV';
% CVtype = '10foldCV';
% CVtype = [CClabel_name 'CCGP'];

WeightSave = 0;
% WeightSave = 1;

% LOOV = 1;
% LOOV = 0;

alpha = 0.0001;


% ========================= save setting =========================
% save_path = 'C:\Users\User\Desktop\201020test\fitclinear_shattering';

% save_path = 'C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results';
% save_path = 'C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results\test';
save_path = 'C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results\fitclinear_shattering';
% save_path = 'C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results\N_vox effect verification';

% save_path = '//143.248.30.72/Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\temp_results\fitclinear_shattering';

% save_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research/Dim_control_PFC_metaRL/results/fitclinear_shattering';
% save_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research/Dim_control_PFC_metaRL/results/fitclinear_shattering_vox_us';
% save_path = '/home/ydsung/A_Research/Dim_control_PFC_metaRL/results/fitclinear_shattering';
% save_path = '/home/ydsung/A_Research/Dim_control_PFC_metaRL/results/fitclinear_shattering_shuffle';

if WeightSave
    basic_name = 'fitclinear_weights';
else
    basic_name = 'fitclinear_accbox';
end
if ~isempty(CClabel_name)
    save_name = [basic_name '_' CClabel_name 'CCGP' exp_sfx];
else
    save_name = [basic_name exp_sfx];
end
% save_name = ['fitclinear_CCGP' exp_sfx CClabel_name];
% save_name = 'fitclinear_accbox_percent';
% save_name = 'fitclinear_accbox_strict';
% save_name = 'fitclinear_accbox_Heo';
% save_name = 'fitclinear_accbox_Kim';
% save_name = 'fitclinear_accbox_Kim_wa';
% save_name = 'fitclinear_accbox_Shin_swar';



% main ====================================================================

fprintf('%d ================================================================================================ \n', id)

% ===== loop1: label load =====
labels = [];

for labi = 1:length(LABEL)
    
    var_name = erase(LABEL{labi}, ' shuff');
    var_name = erase(var_name, ' LOROV');
    var_name = erase(var_name, ' 5fold');
    
    if contains(var_name, ' run')
        runs_temp = var_name((strfind(var_name, ' run')+3):end);
        runs = nan(1,length(runs_temp));
        for i=1:length(runs); runs(i) = str2double(runs_temp(i)); end
        session = arbMBMF_load_var(exp, 'Session', id, []);
        selector = 1 * ismember(session, runs);
        var_name = erase(var_name, var_name(strfind(var_name, ' run'):end));
    else
        selector = [];
    end
    
    % labelling by task variable
    switch var_name
        case 'GoalCond P1'
            label = arbMBMF_load_var(exp, 'GoalCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=1) = nan; 
        case 'GoalCond P2'
            label = arbMBMF_load_var(exp, 'GoalCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=2) = nan;
        case 'GoalCond P3'
            label = arbMBMF_load_var(exp, 'GoalCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=3) = nan;
        case 'UncCond P1'
            label = arbMBMF_load_var(exp, 'UncCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=1) = nan; 
        case 'UncCond P2'
            label = arbMBMF_load_var(exp, 'UncCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=2) = nan;
        case 'UncCond P3'
            label = arbMBMF_load_var(exp, 'UncCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=3) = nan;
        case 'CmplxCond P1'
            label = arbMBMF_load_var(exp, 'CmplxCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=1) = nan; 
        case 'CmplxCond P2'
            label = arbMBMF_load_var(exp, 'CmplxCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=2) = nan;
        case 'CmplxCond P3'
            label = arbMBMF_load_var(exp, 'CmplxCond', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=3) = nan;
        case 'Interaction'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            switch exp
                case {'Lee2014', 'Heo2018'}
                    label(ismember(label,[1 3])) = 1;
                    label(ismember(label,[2 4])) = 2;
                case 'Kim2019'
                    label(ismember(label,[1 4])) = 1;
                    label(ismember(label,[2 3])) = 2;
            end
        case 'Interaction P1'
            label = arbMBMF_load_var(exp, 'Interaction', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=1) = nan; 
        case 'Interaction P2'
            label = arbMBMF_load_var(exp, 'Interaction', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=2) = nan;
        case 'Interaction P3'
            label = arbMBMF_load_var(exp, 'Interaction', id, []);
            phase = class_series2phase_series(label, 'PhaseNumber', 3);
            label(phase~=3) = nan;
        case 'Unc ChoOpt 0'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);   
            A1 = arbMBMF_load_var(exp, 'A1', id, []);
            ChoOpt_old = arbMBMF_load_var(exp, 'ChoOpt_old', id, []);
            label(ismember(label, [1 2]) | ChoOpt_old == 1 | A1 == 1) = nan;
        case 'Unc ChoOpt 1'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);   
            A1 = arbMBMF_load_var(exp, 'A1', id, []);
            ChoOpt_old = arbMBMF_load_var(exp, 'ChoOpt_old', id, []);
            label(ismember(label, [1 2]) | ChoOpt_old == 0 | A1 == 1) = nan;
        case 'Unc in H'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);   
            label(ismember(label, [1 2])) = nan;
        case 'S3G(6,7,8)'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 9])); error('S3 error'); end
            label(~ismember(blkCond, [1 2])) = nan;
            label(label == 9) = nan;
        case 'RG(10,20,40)'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'R', id, ismember(blkCond, [1 2]));
            % sanity check
            if any(~ismember(unq_elms(label), [0 10 20 40])); error('Reward error'); end
            label(label==0) = nan;
        case 'R(10,20,40)'
            label = arbMBMF_load_var(exp, 'R', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [0 10 20 40])); error('Reward error'); end
            label(label==0) = nan;
        case 'S3 new'
            label = arbMBMF_load_var(exp, 'S3', id, []);
        case 'R'
            label = arbMBMF_load_var(exp, 'R', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [0 10 20 40])); error('Reward error'); end
        case 'Goal(6,7,8) uL'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []);    % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 -1])); error('Goal error'); end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label(label == -1) = nan;
            label(~ismember(uncc, 1)) = nan;
        case 'Goal(6,7,8) uH'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []);    % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 -1])); error('Goal error'); end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label(label == -1) = nan;
            label(~ismember(uncc, 2)) = nan;
        case 'Goal(6,7,8)'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 -1])); error('Goal error'); end
            label(label == -1) = nan;
        case 'Goal(6,7,8) binPMB_G0'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            PMB_G = arbMBMF_load_var(exp, 'PMB_G', id, []);
            binPMB_G = regr2class(PMB_G, 2);
            label(binPMB_G~=0) = nan;
        case 'Goal(6,7,8) binPMB_G1'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            PMB_G = arbMBMF_load_var(exp, 'PMB_G', id, []);
            binPMB_G = regr2class(PMB_G, 2);
            label(binPMB_G~=1) = nan;   
        case 'Goal(6.7.8) binChoOpt_G0'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            gt = arbMBMF_load_var(exp, 'GoalCond', id, []);
            opt_g = arbMBMF_load_var(exp, 'ChoOpt', id, gt==1);
            binopt_G = regr2class(opt_g, 2);
            label(binopt_G~=0) = nan;
        case 'Goal(6.7.8) binChoOpt_G1'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            gt = arbMBMF_load_var(exp, 'GoalCond', id, []);
            opt_g = arbMBMF_load_var(exp, 'ChoOpt', id, gt==1);
            binopt_G = regr2class(opt_g, 2);
            label(binopt_G~=1) = nan;
        case 'Goal(6,8) binChoOpt_G0'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            gt = arbMBMF_load_var(exp, 'GoalCond', id, []);
            opt_g = arbMBMF_load_var(exp, 'ChoOpt', id, gt==1);
            binopt_G = regr2class(opt_g, 2);
            label(binopt_G~=0) = nan; label(label==7) = nan;
        case 'Goal(6,8) binChoOpt_G1'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            gt = arbMBMF_load_var(exp, 'GoalCond', id, []);
            opt_g = arbMBMF_load_var(exp, 'ChoOpt', id, gt==1);
            binopt_G = regr2class(opt_g, 2);
            label(binopt_G~=1) = nan; label(label==7) = nan;
        case 'Goal uL'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []);    % 1:low, 2:high
            % sanity check
            switch exp
                case {'Lee2014', 'Heo2018'}
                    if any(~ismember(unq_elms(label), [6 7 8 -1])); error('Goal error'); end
                case 'Kim2019'
                    if any(~ismember(unq_elms(label), [1 2 3])); error('Goal error'); end
            end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label(~ismember(uncc, 1)) = nan;
        case 'Goal uH'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []);    % 1:low, 2:high
            % sanity check
            switch exp
                case {'Lee2014', 'Heo2018'}
                    if any(~ismember(unq_elms(label), [6 7 8 -1])); error('Goal error'); end
                case 'Kim2019'
                    if any(~ismember(unq_elms(label), [1 2 3])); error('Goal error'); end
            end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label(~ismember(uncc, 2)) = nan;
        case 'Goal cL'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            cmplx = arbMBMF_load_var(exp, 'CmplxCond', id, []);    % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(cmplx), [1 2])); error('CmplxCond error'); end
            label(~ismember(cmplx, 1)) = nan;
        case 'Goal cH'
            label = arbMBMF_load_var(exp, 'Goal', id, []);
            cmplx = arbMBMF_load_var(exp, 'CmplxCond', id, []);    % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(cmplx), [1 2])); error('CmplxCond error'); end
            label(~ismember(cmplx, 2)) = nan;
        case 'blkCond only goal 6,-1'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            Goal = arbMBMF_load_var(exp, 'Goal', id, []);
            label(~ismember(Goal, [6 -1])) = nan;
        case 'blkCond in S3 7'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            S3 = arbMBMF_load_var(exp, 'S3', id, []);
            label(~ismember(S3, 7)) = nan;
        case 'blkCond in S2 4'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            S2 = arbMBMF_load_var(exp, 'S2', id, []);
            label(~ismember(S2, 4)) = nan;
        case 'blkCond in A2 right'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            A2 = arbMBMF_load_var(exp, 'A2', id, []);
            label(~ismember(A2, 2)) = nan;
        case 'Goal(1,2,3) x Unc' % Kim2019 only
            goal = arbMBMF_load_var(exp, 'Goal', id, []);    % 1 2 3
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []); % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(goal), [1 2 3])); error('Goal error'); end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label = NaN(1,length(goal)); pairs = [goal; uncc];
            classmap = containers.Map({'1  1', '2  1', '3  1', ...
                                       '1  2', '2  2', '3  2'}, ...
                                       1:6);
            for t=1:length(goal)
                if ~isnan(goal(t))
                    label(t) = classmap(num2str(pairs(:,t)'));
                end
            end
        case 'Goal(6,7,8) x Unc'
            goal = arbMBMF_load_var(exp, 'Goal', id, []);    % 6 7 8 -1
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []); % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(goal), [6 7 8 -1])); error('Goal error'); end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label = NaN(1,length(goal)); pairs = [goal; uncc];
            classmap = containers.Map({'6  1', '7  1', '8  1', ...
                                       '6  2', '7  2', '8  2'}, ...
                                       1:6);
            for t=1:length(goal)
                if ~any(ismember(pairs(:,t)', -1))
                    label(t) = classmap(num2str(pairs(:,t)'));
                end
            end
        case 'Goal x Unc'
            goal = arbMBMF_load_var(exp, 'Goal', id, []);    % 6 7 8 -1
            uncc = arbMBMF_load_var(exp, 'UncCond', id, []); % 1:low, 2:high
            % sanity check
            if any(~ismember(unq_elms(goal), [6 7 8 -1])); error('Goal error'); end
            if any(~ismember(unq_elms(uncc), [1 2])); error('UncCond error'); end
            label = NaN(1,length(goal)); pairs = [goal; uncc];
            classmap = containers.Map({'6  1', '7  1', '8  1', '-1  1', ...
                                       '6  2', '7  2', '8  2', '-1  2'}, ...
                                       1:8);
            for t=1:length(goal)
                label(t) = classmap(num2str(pairs(:,t)'));
            end
        case 'blkCond binChoOpt2 perC 0'
            var2 = arbMBMF_load_var(exp, 'blkCond', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2); if any(~ismember(var2set, [1 2 3 4])); error('err'); end
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan;
        case 'blkCond binChoOpt2 perC 1'
            var2 = arbMBMF_load_var(exp, 'blkCond', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2); if any(~ismember(var2set, [1 2 3 4])); error('err'); end
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan;
        case 'S3(6,7,9) binChoOpt2 perC 0'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan; label(var2==8) = nan; 
        case 'S3(6,7,9) binChoOpt2 perC 1'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan; label(var2==8) = nan; 
        case 'S3 binChoOpt2 perC 0'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan;
        case 'S3 binChoOpt2 perC 1'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []); ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binChoOpt_temp = regr2class(ChoOpt2(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan;
        case 'blkCond ChoCon 0'
            strChoCon = arbMBMF_load_var(exp, 'StrictChoCon', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, []); label(strChoCon==1) = nan;
        case 'blkCond ChoCon 1'
            strChoCon = arbMBMF_load_var(exp, 'StrictChoCon', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, []); label(strChoCon==0) = nan;       
        case 'currBlk blkchng'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; currBlk = nonzeros(Blk)';
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            % sanity check
            if length(blkLens)~=length(currBlk); error('currBlk problem'); end
            for i = 1:length(blkLens)
                label(chng_idx(i)) = currBlk(i);
            end
        case 'currBlk blkchng+1T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; currBlk = nonzeros(Blk)';
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 2
                    label(chng_idx(i)+1) = currBlk(i);
                end
            end
        case 'currBlk blkchng+2T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; currBlk = nonzeros(Blk)';
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 3
                    label(chng_idx(i)+2) = currBlk(i);
                end
            end
        case 'currBlk blkchng+3T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; currBlk = nonzeros(Blk)';
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 4
                    label(chng_idx(i)+3) = currBlk(i);
                end
            end
        case 'currBlk blkchng+4T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; currBlk = nonzeros(Blk)';
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 5
                    label(chng_idx(i)+4) = currBlk(i);
                end
            end
        case 'prevBlk blkchng'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; prevBlk = [nan nonzeros(Blk)'];
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            % sanity check
            if length(blkLens)+1~=length(prevBlk); error('prevBlk problem'); end
            for i = 1:length(blkLens)
                label(chng_idx(i)) = prevBlk(i);
            end
        case 'prevBlk blkchng+1T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; prevBlk = [nan nonzeros(Blk)'];
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 2
                    label(chng_idx(i)+1) = prevBlk(i);
                end
            end
        case 'prevBlk blkchng+2T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; prevBlk = [nan nonzeros(Blk)'];
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 3
                    label(chng_idx(i)+2) = prevBlk(i);
                end
            end
        case 'prevBlk blkchng+3T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; prevBlk = [nan nonzeros(Blk)'];
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 4
                    label(chng_idx(i)+3) = prevBlk(i);
                end
            end
        case 'prevBlk blkchng+4T'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = nan * ones(size(blkCond));
            blkChange = [1 diff(blkCond)~=0];
            Blk = blkCond .* blkChange; prevBlk = [nan nonzeros(Blk)'];
            chng_idx = [find(blkChange) length(blkCond)];
            blkLens = diff(chng_idx);
            for i = 1:length(blkLens)
                if blkLens(i) >= 5
                    label(chng_idx(i)+4) = prevBlk(i);
                end
            end
        case 'blkCond prevT'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            label(2:end) = label(1:end-1); label(1) = nan;
        case 'blkCond nextT'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            label(1:end-1) = label(2:end); label(end) = nan;
        case 'blkCond_sess1'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, session==1);
        case 'blkCond_sess2'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, session==2);
        case 'blkCond_sess3'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, session==3);
        case 'blkCond_sess4'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, session==4);
        case 'blkCond_sess5'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            label = arbMBMF_load_var(exp, 'blkCond', id, session==5);
            %         case 'GoalCond'
            %             blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            %             label = ismember(blkCond, [1 2]);
        case 'S3_blk1'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, blkCond==1);
        case 'S3_blk2'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, blkCond==2);
        case 'binRelMB'
            label = regr2class(arbMBMF_load_var(exp, 'relMB', id, []), 2);
        case 'binRelMF'
            label = regr2class(arbMBMF_load_var(exp, 'relMF', id, []), 2);
        case 'binRelMAX'
            label = regr2class(arbMBMF_load_var(exp, 'relMAX', id, []), 2);            
        case 'binrelMB3xbinrelMF3'
            relMB3 = arbMBMF_load_var(exp, 'relMB3', id, []);
            binrelMB3 = regr2class(relMB3, 2);
            relMF3 = arbMBMF_load_var(exp, 'relMF3', id, []);
            binrelMF3 = regr2class(relMF3, 2);
            label = bin_labels_combination([binrelMB3; binrelMF3]);
        case 'binrelMAX3xrelComp3'
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, []);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, []);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3_sess1'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, session==1);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, session==1);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3_sess2'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, session==2);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, session==2);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3_sess3'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, session==3);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, session==3);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3_sess4'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, session==4);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, session==4);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3_sess5'
            session = arbMBMF_load_var(exp, 'Session', id, []);
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, session==5);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, session==5);
            label = bin_labels_combination([binrelMAX3; relComp3]);
        case 'binrelMAX3xrelComp3 early'
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, []);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, []);
            label = bin_labels_combination([binrelMAX3; relComp3]);
            label(ceil(length(label)/2):end) = nan;
        case 'binrelMAX3xrelComp3 late'
            relMAX3 = arbMBMF_load_var(exp, 'relMAX3', id, []);
            binrelMAX3 = regr2class(relMAX3, 2);
            relComp3 = arbMBMF_load_var(exp, 'relComp3', id, []);
            label = bin_labels_combination([binrelMAX3; relComp3]);
            label(1:floor(length(label)/2)) = nan;
        case 'blkCond_binChoOpt0'
            ChoOpt = arbMBMF_load_var(exp, 'ChoOpt', id, []);
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            if length(ChoOpt)~=length(blkCond); error('error'); end
            ChoOpt_blk1 = ChoOpt(ismember(blkCond, 1));
            ChoOpt_blk2 = ChoOpt(ismember(blkCond, 2));
            ChoOpt_blk3 = ChoOpt(ismember(blkCond, 3));
            ChoOpt_blk4 = ChoOpt(ismember(blkCond, 4));
            binChoOpt_blk1 = regr2class(ChoOpt_blk1, 2);
            binChoOpt_blk2 = regr2class(ChoOpt_blk2, 2);
            binChoOpt_blk3 = regr2class(ChoOpt_blk3, 2);
            binChoOpt_blk4 = regr2class(ChoOpt_blk4, 2);
            binChoOpt_perblkCond = nan * ones(size(ChoOpt));
            binChoOpt_perblkCond(ismember(blkCond, 1)) = binChoOpt_blk1;
            binChoOpt_perblkCond(ismember(blkCond, 2)) = binChoOpt_blk2;
            binChoOpt_perblkCond(ismember(blkCond, 3)) = binChoOpt_blk3;
            binChoOpt_perblkCond(ismember(blkCond, 4)) = binChoOpt_blk4;
            label = blkCond;
            label(binChoOpt_perblkCond==1) = nan;
        case 'blkCond_binChoOpt1'
            ChoOpt = arbMBMF_load_var(exp, 'ChoOpt', id, []);
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            if length(ChoOpt)~=length(blkCond); error('error'); end
            ChoOpt_blk1 = ChoOpt(ismember(blkCond, 1));
            ChoOpt_blk2 = ChoOpt(ismember(blkCond, 2));
            ChoOpt_blk3 = ChoOpt(ismember(blkCond, 3));
            ChoOpt_blk4 = ChoOpt(ismember(blkCond, 4));
            binChoOpt_blk1 = regr2class(ChoOpt_blk1, 2);
            binChoOpt_blk2 = regr2class(ChoOpt_blk2, 2);
            binChoOpt_blk3 = regr2class(ChoOpt_blk3, 2);
            binChoOpt_blk4 = regr2class(ChoOpt_blk4, 2);
            binChoOpt_perblkCond = nan * ones(size(ChoOpt));
            binChoOpt_perblkCond(ismember(blkCond, 1)) = binChoOpt_blk1;
            binChoOpt_perblkCond(ismember(blkCond, 2)) = binChoOpt_blk2;
            binChoOpt_perblkCond(ismember(blkCond, 3)) = binChoOpt_blk3;
            binChoOpt_perblkCond(ismember(blkCond, 4)) = binChoOpt_blk4;
            label = blkCond;
            label(binChoOpt_perblkCond==0) = nan;
        case 'S3(6,7,9)_binChoOpt2_0'
            ChoOpt1n2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []);
            binChoOpt2 = regr2class(ChoOpt1n2(2,:), 2);
            S3_679 = arbMBMF_load_var(exp, 'S3', id, []);
            S3_679(S3_679==8) = nan;
            label = S3_679; label(binChoOpt2==1) = nan;
        case 'S3(6,7,9)_binChoOpt2_1'
            ChoOpt1n2 = arbMBMF_load_var(exp, 'ChoOpt1n2', id, []);
            binChoOpt2 = regr2class(ChoOpt1n2(2,:), 2);
            S3_679 = arbMBMF_load_var(exp, 'S3', id, []);
            S3_679(S3_679==8) = nan;
            label = S3_679; label(binChoOpt2==0) = nan;
        case 'S3(6,7,9)_binChoOpt_0'
            ChoOpt = arbMBMF_load_var(exp, 'ChoOpt', id, []);
            binChoOpt = regr2class(ChoOpt, 2);
            S3_679 = arbMBMF_load_var(exp, 'S3', id, []);
            S3_679(S3_679==8) = nan;
            label = S3_679; label(binChoOpt==1) = nan;
        case 'S3(6,7,9)_binChoOpt_1'
            ChoOpt = arbMBMF_load_var(exp, 'ChoOpt', id, []);
            binChoOpt = regr2class(ChoOpt, 2);
            S3_679 = arbMBMF_load_var(exp, 'S3', id, []);
            S3_679(S3_679==8) = nan;
            label = S3_679; label(binChoOpt==0) = nan;
        case 'blkCond_binrelMAX_0'
            relMAX = arbMBMF_load_var(exp, 'relMAX', id, []);
            binrelMAX = regr2class(relMAX, 2);
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = blkCond; label(binrelMAX==1) = nan;
        case 'blkCond_binrelMAX_1'
            relMAX = arbMBMF_load_var(exp, 'relMAX', id, []);
            binrelMAX = regr2class(relMAX, 2);
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = blkCond; label(binrelMAX==0) = nan;
        case 'S3(6,7,9)'
            label = arbMBMF_load_var(exp, 'S3', id, []);
            label(label==8) = nan;
        case 'S3(6,7,9)_G'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [1 2]));
            label(label==8) = nan;
        case 'S3(6,7,9)_H'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [3 4]));
            label(label==8) = nan;
        case 'S3(6,7,9)_G new'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [1 2]));
            label(label==8) = nan;
        case 'S3(6,7,9)_H new'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [3 4]));
            label(label==8) = nan;
        case 'S(6,7,9)_G uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, 2));
            label(label==8) = nan;
        case 'S(6,7,9)_H uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, 3));
            label(label==8) = nan;
        case 'R(0,20,40)_G'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'R', id, ismember(blkCond, [1 2]));
            label(label==10) = nan;
        case 'R(0,20,40)_H'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'R', id, ismember(blkCond, [3 4]));
            label(label==10) = nan;
        case 'R(0,20,40)_G uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'R', id, ismember(blkCond, 2));
            label(label==10) = nan;
        case 'R(0,20,40)_H uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []); 
            label = arbMBMF_load_var(exp, 'R', id, ismember(blkCond, 3));
            label(label==10) = nan;
        case 'S3(6,7,9)_uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [2 3]));
            label(label==8) = nan;
        case 'S3(6,7,8)_uL'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, ismember(blkCond, [1 4]));
            label(label==9) = nan;
        case 'S3G(6,7,8)_uL'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, blkCond==1);
            label(label==9) = nan;
        case 'S3G(6,7,8)_uH'
            blkCond = arbMBMF_load_var(exp, 'blkCond', id, []);
            label = arbMBMF_load_var(exp, 'S3', id, blkCond==2);
            label(label==9) = nan;
        case 'S3 early'
            label = arbMBMF_load_var(exp, 'S3', id, []);
            label(ceil(length(label)/2):end) = nan;
        case 'S3 late'
            label = arbMBMF_load_var(exp, 'S3', id, []);
            label(1:floor(length(label)/2)) = nan;
        case 'blkCond early'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            label(ceil(length(label)/2):end) = nan;
        case 'blkCond late'
            label = arbMBMF_load_var(exp, 'blkCond', id, []);
            label(1:floor(length(label)/2)) = nan;
        case 'S3(6,7,9) binPMB 0'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'S3', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 9])); error('S3 error'); end
            label(label == 8) = nan; label(binPMB==1) = nan;
        case 'S3(6,7,9) binPMB 1'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'S3', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [6 7 8 9])); error('S3 error'); end
            label(label == 8) = nan; label(binPMB==0) = nan;
        case 'S3(6,7,9) binPMB perC 0'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan; label(var2==8) = nan; 
        case 'S3(6,7,9) binPMB perC 1'
            var2 = arbMBMF_load_var(exp, 'S3', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan; label(var2==8) = nan; 
        case 'R binPMB perC 0'
            var2 = arbMBMF_load_var(exp, 'R', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan; 
        case 'R binPMB perC 1'
            var2 = arbMBMF_load_var(exp, 'R', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan; 
        case 'R(0,20,40) binPMB 0'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'R', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [0 10 20 40])); error('R error'); end
            label(label == 10) = nan; label(binPMB==1) = nan;
        case 'R(0,20,40) binPMB 1'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'R', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), [0 10 20 40])); error('R error'); end
            label(label == 10) = nan; label(binPMB==0) = nan;
        case 'R(0,20,40) binPMB perC 0'
            var2 = arbMBMF_load_var(exp, 'R', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==1) = nan; label(var2==10) = nan; 
        case 'R(0,20,40) binPMB perC 1'
            var2 = arbMBMF_load_var(exp, 'R', id, []); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            PMB = arbMBMF_load_var(exp, 'PMB', id, []); 
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                binPMB_temp = regr2class(PMB(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binPMB_temp==1);
                var1(temp_opt_idx) = 1;
            end
            label = var2; label(var1==0) = nan; label(var2==10) = nan; 
        case 'Coin(1,2,3) binPMB 0'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'Coin', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), 1:4)); error('Coin error'); end
            label(label == 4) = nan; label(binPMB==1) = nan;
        case 'Coin(1,2,3) binPMB 1'
            PMB = arbMBMF_load_var(exp, 'PMB', id, []);
            binPMB = regr2class(PMB, 2);
            label = arbMBMF_load_var(exp, 'Coin', id, []);
            % sanity check
            if any(~ismember(unq_elms(label), 1:4)); error('Coin error'); end
            label(label == 4) = nan; label(binPMB==0) = nan;
        case 'S3(6,7,9) binQarb2 0'
            label = arbMBMF_load_var(exp, 'S3', id, []);
            Qarb2 = arbMBMF_load_var(exp, 'Qarb2', id, []);
            binQarb2 = regr2class(Qarb2, 2);
            label(label == 8) = nan; label(binQarb2 == 1) = nan;
        case 'S3(6,7,9) binQarb2 1'
            label = arbMBMF_load_var(exp, 'S3', id, []);
            Qarb2 = arbMBMF_load_var(exp, 'Qarb2', id, []);
            binQarb2 = regr2class(Qarb2, 2);
            label(label == 8) = nan; label(binQarb2 == 0) = nan;            
        otherwise
            label = arbMBMF_load_var(exp, var_name, id, selector);
    end
        
    if lab_shuffle || contains(LABEL{labi}, 'shuff')
        valid_lab_idx = ~isnan(label);
        valid_lab_shuffled = label(valid_lab_idx);
        valid_lab_shuffled = valid_lab_shuffled(randperm(length(valid_lab_shuffled)));
        label(valid_lab_idx) = valid_lab_shuffled;
    end
    
    if ~isempty(categorize{labi})
        label = regr2class(label, categorize{labi});
    end
    
    labels = [labels; label];
    
end % labi

% For CCGP: cross-condition labeling for each trial
if ~isempty(CClabel_name)
    CClabel = arbMBMF_load_var(exp, CClabel_name, id, []);
else 
    CClabel = [];
end

% ===== loop2: decoding label for each ROI & event =====

for roii = 1:length(ROI)
    roi = ROI{roii};
    fprintf([roi ' ============================================= \n'])
    
    if ~isempty(ROI_types)
        roi_type = ROI_types{roii};
    end
    
    if contains(roi, 'us')  % for voxel undersampling analysis
        k = strfind(roi, 'us');
        % roi(1:k-1): roi_name
        % roi(k+2): underset id
        % sub_vox: voxel index subsets
        MVPfull = arbMBMF_boldpat(exp, roi(1:k-1), id, ... % 3D (N_voxel x N_event x N_trial)
            'SignalType', SignalType, ...
            'roi_type', roi_type, ...
            'Event', EVENT, ...
            'strict', strict, ...
            'interpolation', interpolation);
        MVPfull = MVPfull(sub_vox{str2double(roi(k+2:end))},:,:);
%         disp(sub_vox{str2double(roi(k+2))})
        fprintf('%d ', sub_vox{str2double(roi(k+2:end))})
        disp(' ')
    else
        switch roi
            case 'ilPFC'
                % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
                lilPFC = arbMBMF_boldpat(exp, 'lilPFC', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                rilPFC = arbMBMF_boldpat(exp, 'rilPFC', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                MVPfull = cat(1, lilPFC, rilPFC);
                
            case 'V1'
                lV1 = arbMBMF_boldpat(exp, 'lV1', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                rV1 = arbMBMF_boldpat(exp, 'rV1', id, ...
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation);
                MVPfull = cat(1, lV1, rV1);
                
            otherwise
                % multi-voxel pattern MVP: [N_voxel x N_event x N_trial]
                MVPfull = arbMBMF_boldpat(exp, roi, id, ... 
                    'SignalType', SignalType, ...
                    'roi_type', roi_type, ...
                    'MaskNum', MaskNum_map(roi), ...
                    'Event', EVENT, ...
                    'strict', strict, ...
                    'interpolation', interpolation, ...
                    'prefix', prefix);
        end
    end

    
    for labi = 1:length(LABEL)
        
        fprintf([LABEL{labi} ' ================ \n'])
        
        % loaded label
        label = labels(labi, :);
        
        for evi = 1:length(EVENT)
            
            % event specific MVP
            if ndims(MVPfull) > 2
                MVP = squeeze(MVPfull(:,evi,:));
            else
                MVP = MVPfull;
            end
            % event name for saving
            if iscell(EVENT)
                ev_name = EVENT{evi};
            else
                ev_name = EVENT(evi);
            end
            
            % sanity check
            trial_len = size(label, 2);
            if size(MVP,2)~=trial_len
                disp('label-pattern size unmatch'); 
                continue
            end
            unqSet = unq_elms(label);
            if length(unqSet)~=N_allclass(labi) 
                disp('label element error'); 
                continue
            end
            
            % linear shattering
            fprintf('%d ', id)
            if isempty(KFold)
                acc_box = linear_shattering(MVP, label, ...
                    'Exp', exp, ...
                    'ID', id, ...
                    'CVtype', CVtype, ...
                    'CClabel', CClabel, ...
                    'Repeat', Repeat, ...
                    'rRepeat', rRepeat, ...
                    'WeightSave', WeightSave, ...
                    'Alpha', alpha, ...
                    'AllClassOnly', 1, ...
                    'verbose', 1);
            else
                acc_box = linear_shattering(MVP, label, ...
                    'Kfold', KFold, ...
                    'CClabel', CClabel, ...
                    'Repeat', Repeat, ...
                    'rRepeat', rRepeat, ...
                    'WeightSave', WeightSave, ...
                    'Alpha', alpha, ...
                    'AllClassOnly', 1, ...
                    'verbose', 1);
            end
            
            % running information
            run_info.exp = exp;
            run_info.ROI = ROI;
            run_info.LABEL = LABEL;
            run_info.N_allclass = N_allclass;
            run_info.categorize = categorize;
            run_info.lab_shuffle = lab_shuffle;
            run_info.MVP_setting.SignalType = SignalType;
            run_info.MVP_setting.roi_type = roi_type;
            run_info.MVP_setting.EVENT = EVENT;
            run_info.MVP_setting.strict = strict;
            run_info.MVP_setting.interpolation = interpolation;
            run_info.MVP_setting.prefix = prefix;
            run_info.stat_setting.Repeat = Repeat;
            run_info.stat_setting.rRepeat = rRepeat;
            run_info.stat_setting.CVtype = CVtype;
            run_info.stat_setting.KFold = KFold;
            run_info.stat_setting.WeightSave = WeightSave;
            run_info.stat_setting.Alpha = alpha;
%             run_info.fitcecoc_setting.Coding = Coding;
%             run_info.fitcecoc_setting.Learners = Learners;
%             run_info.fitcecoc_setting.LOOV = LOOV;
                        
            % save information
            save_info.id = id;
            save_info.roi = ROI{roii};
            save_info.roi_type = roi_type;
            save_info.label_name = LABEL{labi};
            save_info.label_true = label;
            save_info.event = ev_name;
            save_info.save_path = save_path;
            save_info.save_name = save_name;
            save_info.acc_box = acc_box;            
            
            % save
            save([save_path '/' roi_type '/' ROI{roii} '/' , ...
                LABEL{labi} '/' num2str(ev_name) '/', ...
                save_name num2str(id) '.mat'], ...
                'run_info', 'save_info')
            disp([roi_type '/' ROI{roii} '/' , ...
                LABEL{labi} '/' num2str(ev_name) '/', ...
                save_name num2str(id) ' saved'])
            
        end % evi
        
    end % labi
    
end % roii

