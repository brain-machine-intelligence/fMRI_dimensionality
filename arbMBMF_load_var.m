
function Variable = arbMBMF_load_var(Exp, var_name, id, tria_cond)

switch Exp
    case 'Lee2014'
        sbjID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};
        % Subject Names
        Sbj_names = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
            'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};  % Lee2014
        try
            base_path = '/home/ydsung/A_Research';
            SBJload = load([base_path, ...
                '/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
            % Session marks
            % sbj_sess_mark : session change mark for MRI slices
%             sess_mark_struct = load([base_path '/Dim_control_PFC_metaRL/session_marks/sess_mark_Lee2014.mat']);
%             sess_mark_Lee2014 = sess_mark_struct.sess_mark_Lee2014;
        catch
            % Loading sbj data
%             base_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/M3_2014subs_ori/SBJ_structure_each_exp_BESTcollectionJan29']);
            SBJload = load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\SBJ_structure_each_exp_BESTcollectionJan29');
            % Session marks
            % sbj_sess_mark : session change mark for MRI slices
%             sess_mark_struct = load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\sess_mark_Lee2014');
%             sess_mark_Lee2014 = sess_mark_struct.sess_mark_Lee2014;
        end
        SBJcell = struct2cell(SBJload);
        SBJs = SBJcell{:};
        % Session : 1~2, 1~3, 1~4 or 1~5
        N_sessions = 5*ones(1,length(Sbj_names)); N_sessions(11)=4; N_sessions(16)=3; N_sessions(20)=2; N_sessions(21)=4;     % Lee2014
        
    case 'Heo2018'
        
%         sbjID = 1:63; sbjID = num2cell(sbjID, 1);
%         SBJ3_load = load('C:\Users\User\Desktop\RL_depression\PLOS\revision\SBJ_structure_1027_for_revision.mat');
        SBJ3_load = load('SBJ_structure_1027_for_revision.mat');
        SBJs = SBJ3_load.SBJ3;
        N_sessions = 4*ones(1,63);      
        
        sbjID = {1, 3, 11, 12, 18, 19, 20, 24, 25, 27, 32, 33, 34, 35, 37, 38, 40, 41, 42, 44, 46, 48, 50, 51, 52, 53, 54, 55};
%         Sbj_names = {'subject001_ksy',[],'subject03_keb','subject004_ksj','subject5_lwc','subject6_sjh','subject7_kch','subject8_kjs', ...
%             'subject9_lks','subject10_ssh','subject11_hjy','subject12_khs','subject13_pjb','subject14_syh','subject15_kik','subject16_jja', ...
%             'subject16_lsh','subj18_kjh','subj19_kny','subj20_jjh','subject21_lsl','subject22_jsh','subject23_yhj','subj24_kjy','subj25_kjs', ...
%             'subject26_cjb','subj27_kkm','subject28_ljh','subject29_cyj','subject30','subj31_ssw','subj32_khj','subj33_ohy','subj34_yhw','subj35_ajs', ...
%             [],'subj37_shi','subj38_lsh',[],'subj40_kkr','subj41_jhj', 'sbj42_ljk',[],'sbj44_ksy', 'sbj45_jej','sbj46_ses','sbj47_hej','sbj48_kdy', ...
%             'subject49_ljy','sbj50_kdh','sbj51_akr', 'sbj52_jkw','sbj53_nyk','sbj54_jyh', 'sbj55_ker'};   % Heo2018
%         try
%             base_path = '/home/ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/SBJ_structure_2014_yd']);
%         catch
%             % Loading sbj data
%             base_path = '\\143.248.30.94\bmlsamba\ydsung/A_Research';
%             SBJload = load([base_path, ...
%                 '/Dim_control_PFC_metaRL/SBJ_structure_2014_yd']);
%         end
%         SBJcell = struct2cell(SBJload);
%         SBJs = SBJcell{:};
%         N_sessions = 4*ones(1,55);      % Heo2018
%         load([base_path '/Dim_control_PFC_metaRL/session_marks/mark_depression.mat']); sbj_sess_mark = mark_depression(sbj_id,:);   % Heo2018

    case 'Kim2019'
        
        sbjID = 1:24;  sbjID([9, 12, 16]) = [];
        sbjID = num2cell(sbjID);
        
        try
            sbj_struct = ... 
                load('C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace/SBJ_structure_Kim2019.mat');
            SBJs = sbj_struct.SBJ;
        catch
            sbj_struct = ... 
                load('/home/bmlshare/complexity/modelRLsource/result_simul/SBJ_structure.mat');
            SBJs = sbj_struct.SBJ;
        end
        
        N_sessions = [5,5,4,6,5,5,6,6,-1,5,5,-1,5,5,5,-1,6,4,5,6,5,5,5,5];
 
end

if id == 0
    
    mean_target_id = nan * ones(1,length(sbjID));
    for idi = 1:length(sbjID)
        temp_tria_cond = [];
        % goalCond
        if tria_cond == 1
            % 200410: temp_tria_cond - for G trials
            temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(temp_tria_cond, [1, 2]);
        end
        if tria_cond == 2
            % 200417: temp_tria_cond - for H trials
            temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(temp_tria_cond, [3, 4]);
        end
        % uncCond
        if tria_cond == 3
            % 200417: temp_tria_cond - for lUnc trials
            temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(temp_tria_cond, [1, 4]);
        end
        if tria_cond == 4
            % 200417: temp_tria_cond - for hUnc trials
            temp_tria_cond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(temp_tria_cond, [2, 3]);
        end
        % blkCond
        if tria_cond == 11
            blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(blkCond, 1);
        end
        if tria_cond == 22
            blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(blkCond, 2);
        end
        if tria_cond == 33
            blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(blkCond, 3);
        end
        if tria_cond == 44
            blkCond = arbMBMF_load_var(Exp, 'blkCond', idi, []);
            temp_tria_cond = ismember(blkCond, 4);
        end
        
        target = arbMBMF_load_var(Exp, var_name, idi, temp_tria_cond);
        if isvector(target)
            mean_target_id(idi) = nanmean(target);
        else
            mean_target_id(idi) = nanmean(nanmean(target));
        end
                
    end
    Variable = mean_target_id;
    
else
    
    switch Exp
        case 'Lee2014'
            SBJ = SBJs{id};
            N_session = N_sessions(sbjID{id});
%             sbj_sess_mark = sess_mark_Lee2014(sbjID{id},:);
            
            % for relMB simplicity
            % relMB
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            relMB(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMB));
            skips = 1;
            relMB(1:skips) = nan;
            relMB = (relMB-nanmin(relMB(skips:end)))/ ...
                (nanmax(relMB(skips:end))-nanmin(relMB(skips:end)));
            relMB = relMB - nanmean(relMB);
            % relMF
            relMF = SBJ.regressor{1,8}.value(8,:);
            relMF(1:3:end) = nan;
%             skips = ceil(0.15 * length(relMF));
            skips = 1;
            relMF(1:skips) = nan;
            relMF = (relMF-nanmin(relMF(skips:end)))/ ...
                (nanmax(relMF(skips:end))-nanmin(relMF(skips:end)));
            relMF = relMF - nanmean(relMF);
            
        case 'Heo2018'
            SBJ = SBJs{sbjID{id}};
            N_session = N_sessions(sbjID{id});
            
        case 'Kim2019'
            SBJ = SBJs{id};
            N_session = N_sessions(sbjID{id});
            
    end
    
    switch var_name
        
        case 0
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = ones(size(reward));

% -------------------------------------------------------------------------
% exp. condition: session
            
        case 'Session' % session
            Variable = [];
            for sess = 1:N_session
                temp_sess = sess * ones(1, size(SBJ.HIST_behavior_info{sess}, 1));
                Variable = [Variable, temp_sess];
            end
            
% -------------------------------------------------------------------------
% exp. condition: stimulus (state)
            
        case 'S2' % 2 3 4 5
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            
        case 'S3' % 6 7 8 9
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            % sanity check
            if ~all(ismember(Variable, [6 7 8 9]))
                error('(arbMBMF_load_var) S3 error');
            end
            
        case 'Coin' % 6 7 8 9
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            switch Exp
                case 'Kim2019'
                    temp = Variable;
                    Variable(ismember(temp, [5 6])) = 1;
                    Variable(ismember(temp, [7 9])) = 2;
                    Variable(ismember(temp, [8 11])) = 3;
                    Variable(ismember(temp, [4 10])) = 4;
            end
            
        case 'S3G'
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==3) = nan;
            Variable(blk_con==4) = nan;
            
        case 10 % S2 in G
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==3) = nan;
            Variable(blk_con==4) = nan;
                    
        case 12 % S2 in H
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,5);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==1) = nan;
            Variable(blk_con==2) = nan;
            
        case 13 % S3 in H
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,6);
                Variable = [Variable, sess_lab'];
            end
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
            end
            Variable(blk_con==1) = nan;
            Variable(blk_con==2) = nan;
            
            
% -------------------------------------------------------------------------
% exp. condition: block condition (context)
            
        case 'blkCond' % block condition 1:G, 2:G', 3:H, 4:H'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = blk_con;
            
        case 'GoalCond' % block condition [1, 2]:G, [3, 4]:H
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = (blk_con > 2.5) + 1; % 1:G, 2:H
            
        case 'UncCond'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    Variable = ismember(blk_con, [2 3]) + 1; % 1:Low, 2:High
                case 'Kim2019'
                    Variable = ismember(blk_con, [3 4]) + 1; % 1:Low, 2:High
            end
            
        case 'CmplxCond'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = ismember(blk_con, [2 4]) + 1; % 1:low, 2:high
            
        case 'Interaction'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = nan(size(blk_con));
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    Variable(ismember(blk_con, [1 3])) = 1;
                    Variable(ismember(blk_con, [2 4])) = 2;
                case 'Kim2019'
                    Variable(ismember(blk_con, [1 4])) = 1;
                    Variable(ismember(blk_con, [2 3])) = 2;
            end
            
        case 'Goal' 
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    Variable = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        Variable = [Variable, sess_lab'];
                    end
            end
            
        case 'Goal(6,7,8)'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    Variable = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        Variable = [Variable, sess_lab'];
                    end
            end
            Variable(Variable == -1) = nan;
            
        % block/goaltype/unc/cmplx/goal switch timing
        case 'blkSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            Variable = [1 diff(blk_con)~=0];
            
        case 'GoaltypeSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = [1 diff(GT)~=0];
            
        case 'UncSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = [1 diff(UC)~=0];
            
        case 'CmplxSwitch'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = [1 diff(CX)~=0];
            
        case 'GoalSwitch'
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    % goal state 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state
                    G = [];
                    for sess = 1:N_session
                        sess_lab = SBJ.HIST_behavior_info{sess}(:,18);
                        G = [G, sess_lab'];
                    end
                case 'Kim2019'
                    % 1: 1st outcome has max value, 2: 2nd, 3: 3rd
                    G = [];
                    for sess = 1:N_session
                        goalVal = SBJ.HIST_behavior_info{sess}(:,18:20);
                        sess_lab = nan(size(goalVal, 1), 1);
                        maxgoal = max(goalVal, [], 2);
                        bin_max = (goalVal == maxgoal);
                        sess_lab(bin_max(:,1)) = 1;
                        sess_lab(bin_max(:,2)) = 2;
                        sess_lab(bin_max(:,3)) = 3;
                        sess_lab(sum(bin_max,2)>1) = nan;
                        G = [G, sess_lab'];
                    end
            end
            Variable = [1 diff(G)~=0];
            
        % phase (e.g. early/late) after block change
        case 'blkPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
%             Variable = class_series2phase_series(blk_con);
            Variable = class_series2phase_series(blk_con, 'PhaseNumber', 3);
        
        case 'GoaltypePhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 3);
            
        case 'GoaltypePhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            
        case 'SpecGoalPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            Variable(GT==0) = nan;
            
        case 'SpecGoalPhase3'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 3);
            Variable(GT==0) = nan;
            
        case 'FlexGoalPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            GT = ismember(blk_con, [1 2]);  % 1: G, 0: H
            Variable = class_series2phase_series(GT, 'PhaseNumber', 2);
            Variable(GT==1) = nan;
            
        case 'UncPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = class_series2phase_series(UC, 'PhaseNumber', 3);
            
        case 'UncPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    UC = ismember(blk_con, [1 4]);  % 1: Low, 0: High
                case 'Kim2019'
                    UC = ismember(blk_con, [1 2]);  % 1: Low, 0: High
            end
            Variable = class_series2phase_series(UC, 'PhaseNumber', 3);
            
        case 'CmplxPhase'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = class_series2phase_series(CX, 'PhaseNumber', 3);
            
        case 'CmplxPhase2'
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            CX = ismember(blk_con, [1 3]);  % 1: Low, 0: High
            Variable = class_series2phase_series(CX, 'PhaseNumber', 3);
            
% -------------------------------------------------------------------------
% exp. condition: episode (S1-A1-S2-A2-S3)
        case 'episode'
            HISTbhv = SBJ.HIST_behavior_info;
            episodes = [];
            for ss = 1:N_session
                episodes = [episodes; num2str(HISTbhv{ss}(:,4:8))];
            end
            
            Variable = nan(1, size(episodes,1));
            epi_set = unique(episodes, 'rows');
            for e = 1:size(epi_set,1)
                Variable(all(episodes==epi_set(e,:),2)) = e;
            end
            
% -------------------------------------------------------------------------
% behavior: action
            
        case 'A1' % A1 1:left, 2:right?
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,7);
                Variable = [Variable, sess_lab'];
            end
            
        case 'A2' % A2 1:left, 2:right?
            Variable = [];
            for sess = 1:N_session
                sess_lab = SBJ.HIST_behavior_info{sess}(:,8);
                Variable = [Variable, sess_lab'];
            end
            
% -------------------------------------------------------------------------
% behavior: performance
            
        case 'Hit'  % Hit
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = zeros(size(reward));
            Variable(reward > 0) = 1;
            
        case 'R'  % Reward
            reward = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                reward = [reward, sess_behav_info(:,16)'];
            end
            Variable = reward;
%             Variable = reward/max(reward);
            
% -------------------------------------------------------------------------
% behavior: choice optimality
            
        case 'ChoOpt_old'  % Choice optimality (Lee2014)
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = is_opt_choice;
            
        case 'ChoOpt' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'likelihood');
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'likelihood');
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
            
            Variable = mean([opt_likelihood1; opt_likelihood2], 1); % [2 x N_trial]
            
        case 'ChoOpt1n2' % choice optimality (Kim2019)
            switch Exp
                case {'Lee2014', 'Heo2018'}
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
                    
                case 'Kim2019'
                    opt_likelihood1 = [];
                    opt_likelihood2 = [];
                    for sess=1:N_session
                        temp_blk_con = SBJ.HIST_block_condition{sess};
                        sess_behav_info = SBJ.HIST_behavior_info{sess};
                        
                        opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                        opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                        for t_sess = 1:size(temp_blk_con, 2)
                            blk = temp_blk_con(2, t_sess);
                            goal = sess_behav_info(t_sess, 18:20);
                            
                            state1 = sess_behav_info(t_sess, 4);
                            action1 = sess_behav_info(t_sess, 7);
                            L1 = opt_normalization_1st2ndstage_complex(state1, action1, goal, blk, 'likelihood');
                            
                            opt_L1_sess(t_sess) = L1;
                            
                            state2 = sess_behav_info(t_sess, 5);
                            action2 = sess_behav_info(t_sess, 8);
                            L2 = opt_normalization_1st2ndstage_complex(state2, action2, goal, blk, 'likelihood');
                            
                            opt_L2_sess(t_sess) = L2;
                        end
                        opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                        opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                        
                        % Sanity check
                        sess_blkcon = sess_behav_info(:,3)';
                        if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                            error('(arbMBMF_load_var) block condition information error');
                        end
                    end
                    Variable = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            end
                        
        case 'ChoOpt1n2G'
            blk_con = [];
            opt_likelihood1 = [];
            opt_likelihood2 = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                
                blk_con = [blk_con, temp_blk_con(2,:)];
                
                opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                for t_sess = 1:size(temp_blk_con, 2)
                    blk = temp_blk_con(2, t_sess);
                    goal = sess_behav_info(t_sess, 18);
                    
                    state1 = sess_behav_info(t_sess, 4);
                    action1 = sess_behav_info(t_sess, 7);
                    L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                    
                    opt_L1_sess(t_sess) = L1;
                    
                    state2 = sess_behav_info(t_sess, 5);
                    action2 = sess_behav_info(t_sess, 8);
                    L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                    
                    opt_L2_sess(t_sess) = L2;
                end
                opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                
                % Sanity check
                sess_blkcon = sess_behav_info(:,3)';
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                    error('(arbMBMF_load_var) block condition information error');
                end
            end
            ChoOpt1n2 = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            G_idx = ismember(blk_con, [1 2]);
            Variable = nan * ones(size(ChoOpt1n2));
            Variable(:,G_idx) = ChoOpt1n2(:,G_idx);
            
        case 'ChoOpt1n2H'
            blk_con = [];
            opt_likelihood1 = [];
            opt_likelihood2 = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                
                blk_con = [blk_con, temp_blk_con(2,:)];
                
                opt_L1_sess = nan * ones(1, size(temp_blk_con, 2));
                opt_L2_sess = nan * ones(1, size(temp_blk_con, 2));
                for t_sess = 1:size(temp_blk_con, 2)
                    blk = temp_blk_con(2, t_sess);
                    goal = sess_behav_info(t_sess, 18);
                    
                    state1 = sess_behav_info(t_sess, 4);
                    action1 = sess_behav_info(t_sess, 7);
                    L1 = opt_normalization_1st2ndstage_exnan(state1, action1, goal, blk);
                    
                    opt_L1_sess(t_sess) = L1;
                    
                    state2 = sess_behav_info(t_sess, 5);
                    action2 = sess_behav_info(t_sess, 8);
                    L2 = opt_normalization_1st2ndstage_exnan(state2, action2, goal, blk);
                    
                    opt_L2_sess(t_sess) = L2;
                end
                opt_likelihood1 = [opt_likelihood1, opt_L1_sess];
                opt_likelihood2 = [opt_likelihood2, opt_L2_sess];
                
                % Sanity check
                sess_blkcon = sess_behav_info(:,3)';
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0
                    error('(arbMBMF_load_var) block condition information error');
                end
            end
            ChoOpt1n2 = [opt_likelihood1; opt_likelihood2]; % [2 x N_trial]
            H_idx = ismember(blk_con, [3 4]);
            Variable = nan * ones(size(ChoOpt1n2));
            Variable(:,H_idx) = ChoOpt1n2(:,H_idx);

% -------------------------------------------------------------------------
% behavior: choice consistency
        case 'ChoConsist1n2'
            
            % ===== state action history =====
            BHVs = [];  % for A1, S2, A2, G save
            for sess = 1:N_session
                BHVmat = SBJ.HIST_behavior_info{sess}';
                % sanity check
                if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                BHVs = [BHVs, BHVmat([7 5 8 18],:)];
            end
            A1 = BHVs(1,:);
            S2 = BHVs(2,:);
            A2 = BHVs(3,:);
            G = BHVs(4,:);
            
%             % ===== when tria_cond given =====
%             if ~isempty(tria_cond)
%                 if length(tria_cond)~=A1; error('tria_cond error'); end
%                 A1 = A1(tria_cond);
%                 S2 = S2(tria_cond);
%                 A2 = A2(tria_cond);
%             end
            
            ChoCon1 = NaN(1,length(A1));
            ChoCon2 = NaN(1,length(A2));
            Goal_set = [6 7 8 -1];
            for gi = 1:length(Goal_set)
                
                g_idx = (G == Goal_set(gi));
                
                % ===== stage 1 choice consistency =====
                A1_temp = A1(g_idx);
                ChoCon1_temp = [0 diff(A1_temp)==0];
                ChoCon1(g_idx) = ChoCon1_temp;
                
                % ===== stage 2 choice consistency =====
                S2_temp = S2(g_idx);
                A2_temp = A2(g_idx);
                ChoCon2_temp = NaN(1,length(A2_temp));
                S2_set = [2 3 4 5];
                for si = 1:length(S2_set)
                    s2_idx = (S2_temp == S2_set(si));
                    if any(s2_idx)
                        a2 = A2_temp(s2_idx);
                        temp_consist2 = [0 diff(a2)==0];
                        ChoCon2_temp(s2_idx) = temp_consist2;
                    end
                end
                ChoCon2(g_idx) = ChoCon2_temp;
                
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoCon1; ChoCon2];
            
        case 'StrictChoCon'
            
            % ===== state action history =====
            A1 = []; S2 = []; A2 = []; G = [];
            for sess = 1:N_session
                BHVmat = SBJ.HIST_behavior_info{sess}';
                
                % sanity check
                if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                
                % S1 is always the same                
                A1 = [A1, BHVmat(7,:)];
                S2 = [S2, BHVmat(5,:)];
                A2 = [A2, BHVmat(8,:)];
                G = [G, BHVmat(18,:)];
            end
            
%             % ===== when tria_cond given =====
%             if ~isempty(tria_cond)
%                 if length(tria_cond)~=A1; error('tria_cond error'); end
%                 A1 = A1(tria_cond);
%                 S2 = S2(tria_cond);
%                 A2 = A2(tria_cond);
%             end
            
            ChoCon1 = NaN(1,length(A1));
            ChoCon2 = NaN(1,length(A2));
            Goal_set = [6 7 8 -1];
            for gi = 1:length(Goal_set)
                
                g_idx = (G == Goal_set(gi));
                
                % ===== stage 1 choice consistency =====
                A1_temp = A1(g_idx);
                ChoCon1_temp = [0 diff(A1_temp)==0];
                ChoCon1(g_idx) = ChoCon1_temp;
                
                % ===== stage 2 choice consistency =====
                S2_temp = S2(g_idx);
                A2_temp = A2(g_idx);
                ChoCon2_temp = NaN(1,length(A2_temp));
                S2_set = [2 3 4 5];
                for si = 1:length(S2_set)
                    s2_idx = (S2_temp == S2_set(si));
                    if any(s2_idx)
                        a2 = A2_temp(s2_idx);
                        temp_consist2 = [0 diff(a2)==0];
                        ChoCon2_temp(s2_idx) = temp_consist2;
                    end
                end
                ChoCon2(g_idx) = ChoCon2_temp;
                
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = all([ChoCon1; ChoCon2], 1) + 0;
            
        case 'ChoCon ignoring context'
            
            % ===== state action history =====
            A1 = []; S2 = []; A2 = []; % G = [];
            for sess = 1:N_session
                BHVmat = SBJ.HIST_behavior_info{sess}';
                
                % sanity check
                if any(~ismember(BHVmat(7,:),[1 2]));       error('A1 error'); end
                if any(~ismember(BHVmat(5,:),[2 3 4 5]));   error('S2 error'); end
                if any(~ismember(BHVmat(8,:),[1 2]));       error('A2 error'); end
                if any(~ismember(BHVmat(18,:),[6 7 8 -1])); error('Goal error'); end
                
                % S1 is always the same
                A1 = [A1, BHVmat(7,:)];
                S2 = [S2, BHVmat(5,:)];
                A2 = [A2, BHVmat(8,:)];
%                 G = [G, BHVmat(18,:)];
            end
            
%             % ===== when tria_cond given =====
%             if ~isempty(tria_cond)
%                 if length(tria_cond)~=A1; error('tria_cond error'); end
%                 A1 = A1(tria_cond);
%                 S2 = S2(tria_cond);
%                 A2 = A2(tria_cond);
%             end
            
            % ===== stage 1 choice consistency =====
            ChoCon1 = [0 diff(A1)==0];
            
            % ===== stage 2 choice consistency =====
            ChoCon2 = NaN(1,length(A2));
            S2_set = [2 3 4 5];
            for si = 1:length(S2_set)
                s2_idx = (S2 == S2_set(si));
                if any(s2_idx)
                    a2_s2 = A2(s2_idx);
                    ChoCon2(s2_idx) = [0 diff(a2_s2)==0];
                end
            end
            
            % sanity check
            if any(isnan(ChoCon2)); error('stage 2 consistency error'); end
            
            % ===== choice consistency with 2 rows =====
            Variable = [ChoCon1; ChoCon2];
            
            
% -------------------------------------------------------------------------
% behavior: response time

        case 'RTA1'
            RTA1 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA1 = [RTA1, sess_behav_info(:,9)'];
            end
            Variable = RTA1;
            
        case 'RTA2'
            RTA2 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA2 = [RTA2, sess_behav_info(:,10)'];
            end
            Variable = RTA2;
            
        case 'RToptA1'
            RTA1 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA1 = [RTA1, sess_behav_info(:,9)'];
            end
            
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = RTA1;
            Variable(~is_opt_choice) = nan;
            
        case 'RToptA2'
            RTA2 = [];
            for sess=1:N_session
                sess_behav_info = SBJ.HIST_behavior_info{sess};
                RTA2 = [RTA2, sess_behav_info(:,10)'];
            end
            
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            Variable = RTA2;
            Variable(~is_opt_choice) = nan;
                        
        case 'EarlyLateDeltaChoOptAfterBlkChange'
            % Len: # of blkChange + 1
            
        case 'TimeUntilFirstOptCho'
            % Len: # of blkChange + 1
            blk_con = []; S2 = []; A2 = []; goal_state = [];
            for sess=1:N_session
                temp_blk_con = SBJ.HIST_block_condition{sess};                 % 2 X session length / 1:blk, 2:blk_con
                blk_con = [blk_con, temp_blk_con(2,:)];         % 1 X total session length (total trial number)
                
                sess_behav_info = SBJ.HIST_behavior_info{sess};             % session length X 18
                sess_blkcon = sess_behav_info(:,3)';
                S2 = [S2, sess_behav_info(:,5)'];
                A2 = [A2, sess_behav_info(:,8)'];
                goal_state = [goal_state, sess_behav_info(:,18)'];
                
                if sum(temp_blk_con(2,:) ~= sess_blkcon) > 0         % Sanity check
                    error('block condition information error'); end
            end
            
            % Finding choice optimality
            is_opt_choice = zeros(size(goal_state));
            
            % specific/low_unc/40goal
            idx_GL40_opt = [find(blk_con==1 & goal_state==6 & S2 == 4 & A2 == 2), ...
                find(blk_con==1 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/low_unc/20goal
            idx_GL20_opt = [find(blk_con==1 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==1 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==1 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/low_unc/10goal
            idx_GL10_opt = [find(blk_con==1 & goal_state==8 & S2 == 2 & A2 == 2), ...
                find(blk_con==1 & goal_state==8 & S2 == 3 & A2 == 1)];
            % specific/high_unc/40goal
            idx_GH40_opt = [find(blk_con==2 & goal_state==6 & S2 == 4 & A2 > 0), ...
                find(blk_con==2 & goal_state==6 & S2 == 5 & A2 == 2)];
            % specific/high_unc/20goal
            idx_GH20_opt = [find(blk_con==2 & goal_state==7 & S2 == 2 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 3 & A2 == 2), ...
                find(blk_con==2 & goal_state==7 & S2 == 4 & A2 == 1), ...
                find(blk_con==2 & goal_state==7 & S2 == 5 & A2 == 1)];
            % specific/high_unc/10goal
            idx_GH10_opt = [find(blk_con==2 & goal_state==8 & S2 == 2 & A2 > 0), ...
                find(blk_con==2 & goal_state==8 & S2 == 3 & A2 == 1)];
            % flexible/high_unc/uni_goal
            idx_FH_opt = [find(blk_con==3 & goal_state==-1 & S2 == 4 & A2 == 1), ...
                find(blk_con==3 & goal_state==-1 & S2 == 5 & A2 == 2)];
            % flexible/low_unc/uni_goal
            idx_FL_opt = [find(blk_con==4 & goal_state==-1 & S2 == 4 & A2 == 2), ...
                find(blk_con==4 & goal_state==-1 & S2 == 5 & A2 == 1)];
            
            is_opt_choice(idx_GL40_opt) = 1; is_opt_choice(idx_GH40_opt) = 1;
            is_opt_choice(idx_GL20_opt) = 1; is_opt_choice(idx_GH20_opt) = 1;
            is_opt_choice(idx_GL10_opt) = 1; is_opt_choice(idx_GH10_opt) = 1;
            is_opt_choice(idx_FH_opt) = 1; is_opt_choice(idx_FL_opt) = 1;
            
            % block change timing
            chng_idx = find([1 diff(blk_con)~=0]~=0);
            T_FirstOpt = nan * ones(size(chng_idx)); % LEN = # of blkChange
            for i = 1:length(chng_idx)
                changed = chng_idx(i);
                opt_after_chng = is_opt_choice(changed:end);
                T_FirstOpt(i) = find(opt_after_chng==1,1,'first');
            end
            Variable = T_FirstOpt;

            
% -------------------------------------------------------------------------
% RL variable: P_MB

        case 'PMB2'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:);
            Variable = PMB(2:3:end);
        case 'PMB3'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:);
            Variable = PMB(3:3:end);
        case 'PMB'  % PMB
            PMB = SBJ.regressor{1,9}.value(7,:);
            Variable = (PMB(1:3:end)+PMB(2:3:end)+PMB(3:3:end))/3;
        case 'PMB_G'
            PMB = SBJ.regressor{1,9}.value(7,:);
            Variable = (PMB(1:3:end)+PMB(2:3:end)+PMB(3:3:end))/3;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end
            
% -------------------------------------------------------------------------
% RL variable: state prediction error

        case 'SPE2' % SPE in stage 2
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(2:3:end);
        case 'SPE3' % SPE in stage 3
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(3:3:end);
        case 'SPE' % average SPE in stage 2 & 3
            SPE = SBJ.regressor{1,1}.value(7,:);
            Variable = SPE(2:3:end) + SPE(3:3:end) / 2;

% -------------------------------------------------------------------------
% RL variable: reward prediction error

        case 'RPE2' % SPE in stage 2
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(2:3:end);
        case 'RPE3' % SPE in stage 3
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(3:3:end);
        case 'RPE' % average SPE in stage 2 & 3
            RPE = SBJ.regressor{1,2}.value(7,:);
            Variable = RPE(2:3:end) + RPE(3:3:end) / 2;

% -------------------------------------------------------------------------
% RL variable: state-action value
            
        case 'Qarb1' % Qarb in stage 1
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(1:3:end);
        case 'Qarb2' % Qarb in stage 2
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(2:3:end);
        case 'Qarb' % average Qarb in stage 1 & 2
            Qarb = SBJ.regressor{1,13}.value(7,:);
            Variable = Qarb(1:3:end) + Qarb(2:3:end) / 2;
            
% -------------------------------------------------------------------------
% RL variable: reliability

        case 'relMB2' % relMB in stage 2
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end);
        case 'relMB3' % relMB in stage 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(3:3:end);
        case 'relMB' % average relMB in stage 2 & 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end) + relMB(3:3:end) / 2;
        case 'relMB_G' % average relMB in stage 2 & 3
            relMB = SBJ.regressor{1,7}.value(8,:)./...
                sum(SBJ.regressor{1,7}.value(7:9,:), 1);
            Variable = relMB(2:3:end) + relMB(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end            
        
        case 'relMF2' % relMF in stage 2
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(2:3:end);
        case 'relMF3' % relMF in stage 3
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(3:3:end);
        case 'relMF' % average relMF in stage 2 & 3
            relMF = SBJ.regressor{1,8}.value(8,:);
            Variable = relMF(2:3:end) + relMF(3:3:end) / 2;
            
        case 'relMAX2' % relMAX in stage 2
            relMAX = nanmax([relMB; relMF], [], 1);
            Variable = relMAX(2:3:end);
        case 'relMAX3' % relMAX in stage 3
            relMAX = nanmax([relMB; relMF], [], 1);
            Variable = relMAX(3:3:end);
        case 'relMAX' % average relMAX in stage 2 & 3
            relMAX = nanmax([relMB; relMF], [], 1);
            Variable = relMAX(2:3:end) + relMAX(3:3:end) / 2;
        case 'relMAX_G'
            relMAX = nanmax([relMB; relMF], [], 1);
            Variable = relMAX(2:3:end) + relMAX(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end
            
        case 'relDiff2'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end);
        case 'relDiff3'
            relDiff = relMB - relMF;
            Variable = relDiff(3:3:end);
        case 'relDiff'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end) + relDiff(3:3:end) / 2;
        case 'relDiff_G'
            relDiff = relMB - relMF;
            Variable = relDiff(2:3:end) + relDiff(3:3:end) / 2;
            blk_con = [];
            for sess=1:N_session
                blk_con = [blk_con, SBJ.HIST_block_condition{sess}(2,:)];
            end
            switch Exp
                case 'Lee2014'
                    Variable(ismember(blk_con, [3 4])) = nan;
            end
            
        case 'relComp2'
            relComp = (relMB > relMF) + 0;
            Variable = relComp(2:3:end);
        case 'relComp3'
            relComp = (relMB > relMF) + 0;
            Variable = relComp(3:3:end);
        case 'relComp'
            relComp = (relMB > relMF) + 0;
            Variable = (relComp(2:3:end) + relComp(3:3:end)) / 2;

% -------------------------------------------------------------------------
% RL variable: learning rate

        case 'lr'
            Variable = SBJ.model_BayesArb.param(8);
                        
    end
    
% -------------------------------------------------------------------------

    if ~isempty(tria_cond)
        if isscalar(tria_cond)
            switch tria_cond
                case 1
                    blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                    Variable(ismember(blkCond, [3, 4])) = nan;
                case 2
                    blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                    Variable(ismember(blkCond, [1, 2])) = nan;
            end
        else
            if isvector(Variable)
                if any(isnan(tria_cond))
                    Variable(isnan(tria_cond)) = nan;
                else
                    Variable(tria_cond==0) = nan;
                end
            else
                if any(isnan(tria_cond))
                    Variable(:, isnan(tria_cond)) = nan;
                else
                    Variable(:, tria_cond==0) = nan;
                end
            end
        end
    end
    
end

