function mri_ids = SARtime_to_MRIidx(sess_HIST_event_info, TR, delay, varargin)
% Convert S, A, R timepoints in 'HIST_event_info' into MRI file (nii) IDs
% Inputs :  sess_HIST_event_info = HIST_event_info{sess} : [8 x 8*N_trial_sess mat], 8 slots for 1 trial (fix1/S1/A1/fix2/S2/A2/fix3/R)
%           TR : MRI repetition time [sec]
%           delay : time delay by hemodynamic response [sec]
% Output : mri_ids [8 x N_trial mat]

% default input parameters
options = struct('Exp', 'Lee2014', ...
                'strict', [], ...
                'forInterpol', [], ...
                'EventPeriod', 8);
% read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(SARtime_to_MRIidx) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(SARtime_to_MRIidx) %s is not a recognized parameter name', pair{1})
    end
end
ep = options.EventPeriod;

switch options.Exp
    case {'Lee2014', 'Heo2018'}
        
        % Input validity
        if size(sess_HIST_event_info, 1) < 8; error('(SARtime_to_MRIidx) input sess_HIST_event_info does not have 8 rows'); end
        if mod(size(sess_HIST_event_info, 2), ep); error('(SARtime_to_MRIidx) input sess_HIST_event_info N_column is not 8*N_trial'); end
        if ~(TR > 0 && delay > 0); error('(SARtime_to_MRIidx) TR and delay must be > 0'); end
        SARs = sess_HIST_event_info(7,:);
        S1s = SARs(2:8:end); if ~all(ismember(S1s,1)); error('(SARtime_to_MRIidx) S1~=1 case exists'); end
        S2s = SARs(5:8:end); if ~all(ismember(S2s,[2 3 4 5])); error('(SARtime_to_MRIidx) S2~=2,3,4,5 case exists'); end
        Rs = SARs(8:8:end); if ~all(ismember(round(Rs),[6 7 8 9])); error('(SARtime_to_MRIidx) R~=6,7,8,9 case exists'); end
        fix1s = SARs(1:8:end); fix2s = SARs(4:8:end); fix3s = SARs(7:8:end);
        if ~all(ismember([fix1s, fix2s, fix3s],0.5)); error('(SARtime_to_MRIidx) fixation~=0.5 case exists'); end
        
        % Event times in session
        tfix1 = sess_HIST_event_info(4,1:8:end);    % first fixation
        tS1 = sess_HIST_event_info(4,2:8:end);      % state 1
        tA1 = sess_HIST_event_info(4,3:8:end);      % action 1
        tfix2 = sess_HIST_event_info(4,4:8:end);    % second fixation
        tS2 = sess_HIST_event_info(4,5:8:end);      % state 2
        tA2 = sess_HIST_event_info(4,6:8:end);      % action 2
        tfix3 = sess_HIST_event_info(4,7:8:end);    % third fixation
        tR = sess_HIST_event_info(4,8:8:end);       % Reward (terminal state)
        
    case 'Kim2019'
        
        % Input validity
        if size(sess_HIST_event_info, 1) < 8; error('(SARtime_to_MRIidx) input sess_HIST_event_info does not have 8 rows'); end
        if ~(TR > 0 && delay > 0); error('(SARtime_to_MRIidx) TR and delay must be > 0'); end
        
        % trial stage (1, 2, 3)
        stage = sess_HIST_event_info(3,:);
        
        % event time series
        SARs = sess_HIST_event_info(7,:);
        % row7 - state.
        % 0.5: fixation mark on,
        % 1~11: S1~S11 >> (+/-)0.1: (with win/lost msg),
        % 21-24:L1/R1/L2/R2,
        % 30: a short blank page display,
        % -99:fail to choose in time limit, (-) when display off
        
        % stage 1 (event 1~3: fix1, S1, A1)
        tfix1 = sess_HIST_event_info(4, ismember(stage, 1) & ...
            ismember(SARs, 0.5));
        tS1 = sess_HIST_event_info(4, ismember(stage, 1) & ...
            ~ismember(SARs, 0.5) & ...
            ismember(round(SARs), 1));
        tA1 = sess_HIST_event_info(4, ismember(stage, 1) & ...
            ismember(round(SARs), [21 22 -99]));
        
        % stage 2 (event 4~6: fix2, S2, A2)
        tfix2 = sess_HIST_event_info(4, ismember(stage, 2) & ...
            ismember(SARs, 0.5));
        tS2 = sess_HIST_event_info(4, ismember(stage, 2) & ...
            ismember(round(SARs), [2 3]));
        tA2 = sess_HIST_event_info(4, ismember(stage, 2) & ...
            ismember(round(SARs), [21:24 -99]));
        
        % stage 3 (event 7~9: fix3, S3, fix3)
        fix3_idx = find(stage==3 & SARs==0.5);
        tfix3 = sess_HIST_event_info(4, fix3_idx(stage(fix3_idx-1)==2));
        tR = sess_HIST_event_info(4, ismember(stage, 3) & ...
            ismember(round(SARs), 4:11));
        tend = sess_HIST_event_info(4, fix3_idx(stage(fix3_idx-1)==3));
end

mri_ids = [tfix1; tS1; tA1; tfix2; tS2; tA2; tfix3; tR];
% ceiling to integer ID (no interpolation) or not (for interpolation)
interpol = options.forInterpol;
if ~isempty(interpol) && interpol~=0 && interpol~=false
    mri_ids = (mri_ids+delay)/TR;
else
    mri_ids = ceil((mri_ids+delay)/TR);
end
% (mainly for no interpolation) strict allocation without MRI vol ID overlap
strict = options.strict;
if ~isempty(strict) && strict~=0 && strict~=false
    next_ids = [mri_ids(2:end,:); [mri_ids(1,2:end), nan]];
    mri_ids(mri_ids == next_ids) = nan;
end

% Output validity
unfold = reshape(mri_ids, 1, []); unfold_idx = ~isnan(unfold);
sorted = sort(unfold); sorted_idx = ~isnan(sorted);
if ~all(unfold(unfold_idx) == sorted(sorted_idx))
    error('(SARtime_to_MRIidx) Output validity problem')
end

end

