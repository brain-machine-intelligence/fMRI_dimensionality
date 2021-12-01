function acc_box = linear_shattering(pattern, label, varargin)
% Neural pattern dimension measurement by implementation of possible SVMs
% (input) pattern : N_vox X N_timepoints
%         label : 1 X N_timepoints (values : 1, 2, ... , c, null label = NaN)
% (output) acc_box
    
%     % Input validity
%     p = inputParser;
%     validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%     addRequired(p,'pattern',@ismatrix);
%     addRequired(p,'label',@isvector);
%     addParameter(p,'Kfold',validScalarPosNum);
%     addParameter(p,'CVtype',validScalarPosNum);
%     addParameter(p,'Repeat',validScalarPosNum);
%     addParameter(p,'rRepeat',validScalarPosNum);
%     addParameter(p,'Alpha',validScalarPosNum);
%     addParameter(p,'AllClassOnly',validScalarPosNum);
%     addParameter(p,'verbose',validScalarPosNum);
%     parse(p,pattern,label,varargin{:});

    % 'label' should be a row vector & pattern size should be matched
    if size(label,1)>1; label=label'; end
    if size(pattern,2)~=size(label,2)
        if size(pattern,1)==size(label,2)
            warning('(liear_shattering) Input pattern transposed');
            pattern = pattern';
        else
            error('(linear_shattering) Neural pattern matrix and label vector must have the same timepoints!')
        end
    end
    
    % Setting default values for options
    options = struct('Exp',[],'ID',[],'Kfold',[],'CVtype',[],'CClabel',[],'Repeat',1,'rRepeat',0, ...
                    'WeightSave', 0, ...
                    'Alpha',0.01,'AllClassOnly',0,'verbose',1);
    optionNames = fieldnames(options);
    
    % Option value matching & name validity
    n_args = length(varargin);
    if round(n_args/2) ~= n_args/2
        fprintf('length(varargin)=%d \n',n_args)
        error('(linear_shattering) varargin pair error')
    end
    for pair = reshape(varargin, 2, [])
        if ismember(pair{1}, optionNames)
            options.(pair{1}) = pair{2};
        else
            error('(linear_shattering) %s is not recognized parameter name', pair{1});
        end
    end
    exp = options.('Exp');
    id = options.('ID');
    Kfold = options.('Kfold');
    CVtype = options.('CVtype');
    CClabel = options.('CClabel');
    Repeat = options.('Repeat');
    rRepeat = options.('rRepeat');
    WeightSave = options.('WeightSave');
    Alpha = options.('Alpha');
    AllClassOnly = options.('AllClassOnly');
    verbose = options.('verbose');
    
    % 'label' must not have 0 as an element
    if ismember(0,label); label=label+1; end
    
    lab_classes = unique(label(~isnan(label)));
    c = length(lab_classes);
    if verbose
        fprintf(['(linear_shattering) Input classes: ' num2str(lab_classes)]);
    end
    
    bin_conds = cell(1, c);
    for i=1:c
        bin_conds{i} = [lab_classes(i) 0];
    end
    pat_size = size(pattern, 1);
    
    acc_box.ShatteredClassesNumber = cell(1, c);
    acc_box.options = options;
    
    % Shattering for 2 ~ c classes
    for m=2:c
        if AllClassOnly
            if m ~= c; continue; end
        end
%         CASE.n_conditions = m;
% %         CASE.Tc = 2^m - 2;
%         CASE.separable = zeros(nchoosek(c,m), 2^m - 2);
%         CASE.separable2 = zeros(nchoosek(c,m), 2^m - 2);
% %         CASE.pca_outputs = cell(nchoosek(c,m), 2^m - 2, Repeat);
%         CASE.PRs = nan * ones(nchoosek(c,m), 2^m - 2, Repeat);
%         CASE.nPCs = nan * ones(nchoosek(c,m), 2^m - 2, Repeat);
% %         CASE.latents = cell(nchoosek(c,m), 2^m - 2);
%         CASE.accs = nan * ones(nchoosek(c,m), 2^m - 2, Repeat);
%         CASE.acc_mean = nan * ones(nchoosek(c,m), 2^m - 2);
%         CASE.acc_ste = nan * ones(nchoosek(c,m), 2^m - 2);
%         CASE.accs_r = nan * ones(nchoosek(c,m), 2^m - 2, Repeat);
%         CASE.acc_mean_r = nan * ones(nchoosek(c,m), 2^m - 2);
%         CASE.acc_ste_r = nan * ones(nchoosek(c,m), 2^m - 2);
%         CASE.positive_labels = cell(nchoosek(c,m), 2^m - 2);
%         CASE.sample_sizes = cell(nchoosek(c,m), 2^m - 2);
%         CASE.posi_or_neg_lab_only = zeros(nchoosek(c,m), 2^m - 2);
%         CASE.insufficient_for_CV = zeros(nchoosek(c,m), 2^m - 2);
%         if rRepeat; CASE.accs_rr = nan * ones(nchoosek(c,m), 2^m - 2, Repeat, rRepeat); end

        CASE.n_conditions = m;
%         CASE.Tc = (2^m - 2)/2;
        CASE.separable = zeros(nchoosek(c,m), (2^m - 2)/2);
        CASE.separable2 = zeros(nchoosek(c,m), (2^m - 2)/2);
%         CASE.pca_outputs = cell(nchoosek(c,m), (2^m - 2)/2, Repeat);
        CASE.PRs = nan * ones(nchoosek(c,m), (2^m - 2)/2, Repeat);
        CASE.nPCs = nan * ones(nchoosek(c,m), (2^m - 2)/2, Repeat);
%         CASE.latents = cell(nchoosek(c,m), (2^m - 2)/2);
        CASE.accs = nan * ones(nchoosek(c,m), (2^m - 2)/2, Repeat);
        CASE.accs_run_mat = cell(nchoosek(c,m), (2^m - 2)/2);
        CASE.acc_mean = nan * ones(nchoosek(c,m), (2^m - 2)/2);
        CASE.acc_ste = nan * ones(nchoosek(c,m), (2^m - 2)/2);
        CASE.accs_r = nan * ones(nchoosek(c,m), (2^m - 2)/2, Repeat);
        CASE.acc_mean_r = nan * ones(nchoosek(c,m), (2^m - 2)/2);
        CASE.acc_ste_r = nan * ones(nchoosek(c,m), (2^m - 2)/2);
        if WeightSave
            CASE.betas = nan(nchoosek(c,m), (2^m - 2)/2, pat_size, Repeat);
            CASE.betas_per_run = cell(nchoosek(c,m), (2^m - 2)/2);
            CASE.beta_mean = nan(nchoosek(c,m), (2^m - 2)/2, pat_size);
        end
        CASE.positive_labels = cell(nchoosek(c,m), (2^m - 2)/2);
        CASE.sample_sizes = cell(nchoosek(c,m), (2^m - 2)/2);
        CASE.posi_or_neg_lab_only = zeros(nchoosek(c,m), (2^m - 2)/2);
        CASE.insufficient_for_CV = zeros(nchoosek(c,m), (2^m - 2)/2);
        if rRepeat; CASE.accs_rr = nan * ones(nchoosek(c,m), (2^m - 2)/2, Repeat, rRepeat); end
                            
        class_combins = nchoosek(1:c, m);       % class_combins : cCm X m
        for com = 1:size(class_combins,1)
            positive_label_set = combvec(bin_conds{class_combins(com,:)});     % positive_label_set : m X 2^m
            p_cnt = 0;
            % for p = 1:size(positive_label_set,2)
            for p = 1:size(positive_label_set,2)/2
                if sum(positive_label_set(:,p)) ~= 0 && length(nonzeros(positive_label_set(:,p))) < m % length(class_combins(com,:))
                    p_cnt = p_cnt+1;
                    new_label = nan*ones(size(label));
                    % positive(1) or negative(0) labelling for chosen classes
                    for selected_cond = lab_classes(class_combins(com, :))
                        new_label(label==selected_cond) = ismember(selected_cond, positive_label_set(:,p));
                    end
                    
%                     disp(positive_label_set(:,p))
                    
                    pos_length = sum(new_label==1); neg_length = sum(new_label==0);
                    if pos_length * neg_length == 0
                        CASE.posi_or_neg_lab_only(com, p_cnt) = 1;
                    else
                        
                        % Comparing repeated undersampled classification result & random shuffled data result
                        accs = nan * ones(1,Repeat);
                        accs_r = nan * ones(1,Repeat);
                        accs_rr = nan * ones(Repeat,rRepeat);
                        
                        if WeightSave; betas = nan(pat_size, Repeat); end
                        
                        if ~isempty(exp)
                            N_sess = length(unq_elms(arbMBMF_load_var(exp, 'Session', id, [])));
                            accs_per_run = nan(Repeat, N_sess);
                            if WeightSave; betas_per_run = nan(Repeat, pat_size, N_sess); end
                        else
                            accs_per_run = [];
                            betas_per_run = [];
                        end
                        for r = 1:Repeat
                            
                            rng('shuffle');
%                             pca_output.coeff = coeff;
%                             pca_output.score = score;
%                             pca_output.latent = latent;
%                             pca_output.explained = explained;
%                             pca_output.PR = sum(latent)^2/sum(latent.^2);
%                             pca_output.nPC = sum( cumsum(latent)/sum(latent) <= 0.9 );
%                             CASE.pca_outputs{com, p_cnt, r} = pca_output;
                            
                            % SVM training
                            if isempty(Kfold)
                                switch CVtype
                                    case 'LOROV'
                                        runs = arbMBMF_load_var(exp, 'Session', id, []);
                                        [runSet, runN] = unq_elms(runs);
                                        accs_run = nan(1, length(runSet));
                                        if WeightSave; betas_run = nan(pat_size, length(runSet)); end
                                        if isempty(CClabel) % dumy label for empty CClabel
                                            temp_CClabel = ones(size(new_label));
                                        else
                                            temp_CClabel = CClabel;
                                        end
                                        [clabSet, clabN] = unq_elms(temp_CClabel);
                                        for run = runSet
                                            % recording CCGP across CClabel conditions
                                            accs_con = nan(1, length(clabSet));
                                            if WeightSave; beta_con = nan(pat_size, length(clabSet)); end
                                            for cci = 1:length(clabSet)
                                                if isempty(CClabel)
                                                    test_idx = (runs == run);
                                                    train_idx = (runs ~= run);
                                                else
                                                    test_idx = (runs == run) & (temp_CClabel == clabSet(cci));
                                                    train_idx = (runs ~= run) & (temp_CClabel ~= clabSet(cci));
                                                end
                                                % Balancing label by undersampling
                                                [label_test, MVP_test] = ...
                                                    lab_pat_undersample(new_label(test_idx), pattern(:, test_idx));
                                                [label_train, MVP_train] = ...
                                                    lab_pat_undersample(new_label(train_idx), pattern(:, train_idx));
                                                Mdl = fitclinear(MVP_train, label_train, ...
                                                    'Learner', 'svm', 'ObservationsIn', 'columns', 'Regularization', 'ridge', 'Solver', 'dual');
                                                accs_con(cci) = 1 - loss(Mdl, MVP_test', label_test');
                                                if WeightSave; beta_con(:, cci) = Mdl.Beta; end
                                            end
                                            % row * col / scalar
                                            accs_run(run) = accs_con * clabN / sum(clabN);
                                            % rows (mat) * col / scalar
                                            if WeightSave; betas_run(:, run) = beta_con * clabN / sum(clabN); end
                                        end
                                        accs(r) = accs_run * runN / sum(runN);
                                        accs_per_run(r,:) = accs_run;
                                        if WeightSave; betas(:, r) = betas_run * runN / sum(runN); end
                                        if WeightSave; betas_per_run(r,:,:) = betas_run; end
                                end
                            else % 'Kfold is not empty'
                                
                                if isempty(CClabel)
                                    
                                    % Balancing label by undersampling
                                    [final_lab, final_pat] = lab_pat_undersample(new_label, pattern);
                                    
                                    
                                    [~,~,latent] = pca(final_pat');
                                    CASE.PRs(com, p_cnt, r) = sum(latent)^2/sum(latent.^2);
                                    CASE.nPCs(com, p_cnt, r) = sum( cumsum(latent)/sum(latent) <= 0.9 );
                                    
                                    if length(final_lab) < 2*Kfold
                                        CASE.insufficient_for_CV(com, p_cnt) = CASE.insufficient_for_CV(com, p_cnt) + 1;
                                    else
                                        rand(floor(mod(sum(clock*10),10000)));
                                        cvsvm = fitclinear(final_pat, final_lab, 'Kfold', Kfold, ...
                                            'Learner', 'svm', 'ObservationsIn', 'columns', 'Regularization', 'ridge', 'Solver', 'dual');
                                        acc = 1 - kfoldLoss(cvsvm);
                                        %                                 prediction = kfoldPredict(cvsvm);       % column vector
                                        %                                 acc = sum(final_lab'==prediction)/length(final_lab);
                                        accs(r) = acc;
                                        
                                        clear cvsvm; clear final_pattern;
                                    end % Kfold
                                    
                                else
                                    % if CClabel is not empty
                                    temp_CClabel = CClabel;
                                    [clabSet, clabN] = unq_elms(temp_CClabel);
                                    if Kfold ~= length(clabSet)
                                        error('(linear shatterin) Kfold and CClabel class # unmatching')
                                    end
                                    
                                    accs_con = nan(1, length(clabSet));
                                    for cci = 1:length(clabSet)
                                        test_idx = (temp_CClabel == clabSet(cci));
                                        train_idx = (temp_CClabel ~= clabSet(cci));
                                        % Balancing label by undersampling
                                        [label_test, MVP_test] = ...
                                            lab_pat_undersample(new_label(test_idx), pattern(:, test_idx));
                                        [label_train, MVP_train] = ...
                                            lab_pat_undersample(new_label(train_idx), pattern(:, train_idx));
                                        rand(floor(mod(sum(clock*10),10000)));
                                        Mdl = fitclinear(MVP_train, label_train, ...
                                            'Learner', 'svm', 'ObservationsIn', 'columns', 'Regularization', 'ridge', 'Solver', 'dual');
                                        accs_con(cci) = 1 - loss(Mdl, MVP_test', label_test');
                                    end
                                    
                                    accs(r) = accs_con * clabN / sum(clabN);
                                    
                                end 
                                
                                if rRepeat
                                    for rr = 1:rRepeat
                                        shuffle_lab = final_lab(randperm(length(final_lab)));
                                        cvsvm_r = fitclinear(final_pat, shuffle_lab, 'Kfold', Kfold, ...
                                            'Learner', 'svm', 'ObservationsIn', 'columns', 'Regularization', 'ridge', 'Solver', 'dual');
                                        acc_r = 1 - kfoldLoss(cvsvm_r);
                                        %                                         prediction_r = kfoldPredict(cvsvm_r);       % column vector
                                        %                                         acc_r = sum(shuffle_lab'==prediction_r)/length(shuffle_lab);
                                        accs_rr(r, rr) = acc_r;
                                        clear cvsvm_r;
                                    end
                                    accs_r(r) = mean(accs_rr(r, :));
                                else
                                    accs_r(r) = 0.5;
                                end
                                
                            end % Kfold empty
                            
                        end % Repeat
                        
                        if ~all(isnan(accs))
                            h = ttest(accs(~isnan(accs)), accs_r(~isnan(accs_r)), 'Tail', 'right', 'Alpha', Alpha);
                            h2 = ttest(accs(~isnan(accs)), 0.5 * ones(1,sum(~isnan(accs))), 'Tail', 'right', 'Alpha', Alpha);
                            if isnan(h); disp(' '); disp([num2str(accs) '/' num2str(accs_r)]); 
                            elseif h; CASE.separable(com,p_cnt)=1; end
                            if h2; CASE.separable2(com,p_cnt)=1; end
                            CASE.accs(com, p_cnt, :) = accs;
                            if ~isempty(exp); CASE.accs_run_mat{com, p_cnt} = accs_per_run; end
%                             if ~isempty(exp); CASE.accs_run_mat(com, p_cnt, :) = accs_run_mat; end
                            CASE.acc_mean(com, p_cnt) = nanmean(accs);
                            CASE.acc_ste(com, p_cnt) = nanstd(accs)/sqrt(sum(~isnan(accs)));
                            CASE.accs_r(com, p_cnt, :) = accs_r;
                            CASE.acc_mean_r(com, p_cnt) = nanmean(accs_r);
                            CASE.acc_ste_r(com, p_cnt) = nanstd(accs_r)/sqrt(sum(~isnan(accs_r)));
                            % if CCGP, beta is not saved
                            if isempty(CClabel)
                                if WeightSave; CASE.betas(com, p_cnt, :, :) = betas; end
                                % memory issue
                                % if ~isempty(exp); CASE.betas_per_run{com, p_cnt} = betas_per_run; end
                                if WeightSave; CASE.beta_mean(com, p_cnt, :) = nanmean(betas, 2); end
                            end
                            if rRepeat; CASE.accs_rr(com, p_cnt, :, :) = accs_rr; end
                        end
                        CASE.insufficient_for_CV(com, p_cnt) = CASE.insufficient_for_CV(com, p_cnt)/Repeat;
                        
                    end % pos_length * neg_length == 0
                    
                    CASE.positive_labels{com, p_cnt} = nonzeros(positive_label_set(:,p));
                    CASE.sample_sizes{com, p_cnt} = [pos_length neg_length];
                end % p_cnt
            end % p
        end % com
        acc_box.ShatteredClassesNumber{m} = CASE;
        if verbose; fprintf(', %d/%d', m, c); end
    end % m
    if verbose; disp(' done'); end
end

