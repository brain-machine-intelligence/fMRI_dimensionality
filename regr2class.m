function [RegrClassVector, separate_cutoff] = regr2class(regressor, sep, varargin)
% Convert continuous regressor vector into class vector
% Inputs :  regressor : regressor time series (row vector)
%           sep : # of class separation
% Output : RegrClassVector : converted class vector
%          separate_cutoff : cutoff values & type [struct]

% Input validity
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'regressor',@isvector);
addRequired(p,'sep',validScalarPosNum);
addParameter(p,'NumberForwardIgnore',validScalarPosNum);
addParameter(p,'NumberBackwardIgnore',validScalarPosNum);
parse(p,regressor,sep,varargin{:});

if mod(length(varargin),2)
    fprintf('length(varargin{:}) = %d', length(length(varargin{:})));
    error('(regr2class) varargin pair error');
end

% Setting default values
options = struct('NumberForwardIgnore',0,'NumberBackwardIgnore',0);

for pair = reshape(varargin,2,[])
    if ismember(pair{1}, fieldnames(options))
        options.(pair{1}) = pair{2};
    else
        error('(regr2class) %s is not recognized parameter name', pair{1});
    end
end

if options.('NumberForwardIgnore'); regressor = regressor(options.('NumberForwardIgnore'):end); end
if options.('NumberBackwardIgnore'); regressor = regressor(1:options.('NumberBackwardIgnore')); end

if iscolumn(regressor)
    regressor = regressor';     % label should be a row vector
end

% regressor into class label
if sep == 2     % binary labelling
    med_boundary = prctile(regressor, 50);
    mean_boundary = nanmean(regressor);
    
    % med_boundary or mean_boundary
    N_1type_med = min(sum(regressor>med_boundary),sum(regressor<med_boundary));
    N_1type_mean = min(sum(regressor>mean_boundary),sum(regressor<mean_boundary));
    if N_1type_med + N_1type_mean == 0
        error('(regr2class) Cutoff error in regressor binary labelling')
    elseif N_1type_med > N_1type_mean
        RegrClassVector = regressor > med_boundary;
        separate_cutoff.type = 'median';
        separate_cutoff.value = med_boundary;
    elseif N_1type_med <= N_1type_mean
        RegrClassVector = regressor > mean_boundary;
        separate_cutoff.type = 'mean';
        separate_cutoff.value = mean_boundary;
    end
    RegrClassVector = double(RegrClassVector);
    RegrClassVector(isnan(regressor)) = nan;
    
else            % multi-level labelling
    
    % percentile-based labelling
    cutoff_bound1 = prctile(regressor, round(100/sep):round(100/sep):100);
    rep_extr_reg1 = repmat(regressor, length(cutoff_bound1)-1, 1);
    reg_label1 = sum(rep_extr_reg1 > cutoff_bound1(1:end-1)', 1) + 1;
    reg_label1(isnan(regressor)) = nan;
    
    % value-based labelling
    cutoff_bound2 = linspace(min(regressor), max(regressor), sep+1);
    rep_extr_reg2 = repmat(regressor, length(cutoff_bound2)-2, 1);
    reg_label2 = sum(rep_extr_reg2 > cutoff_bound2(2:end-1)', 1) + 1;
    reg_label2(isnan(regressor)) = nan;
    
    % Comparing 2 labelling & choosing richer one
    [~, ~, ic1] = unique(reg_label1(~isnan(reg_label1)));
    uni_counts1 = accumarray(ic1, 1); % counting unique values
    [~, ~, ic2] = unique(reg_label2(~isnan(reg_label2)));
    uni_counts2 = accumarray(ic2, 1); % counting unique values
    
    if min(uni_counts1) + min(uni_counts2) == 0
        error('(regr2class) Cutoff error in regressor multi-level labelling')
    elseif min(uni_counts1) > min(uni_counts2)
        RegrClassVector = reg_label1;
        separate_cutoff.type = ['percentile' num2str(sep)];
        separate_cutoff.value = cutoff_bound1;
    elseif min(uni_counts1) <= min(uni_counts2)
        RegrClassVector = reg_label2;
        separate_cutoff.type = ['linspace' num2str(sep)];
        separate_cutoff.value = cutoff_bound2;
    end
    
end

