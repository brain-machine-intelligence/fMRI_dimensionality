function [new_label, new_pattern, N_sampling] = lab_pat_undersample(label, pattern, varargin)
%Balance labels and corresponding patterns(features) by undersampling given
%label vector and pattern matrix. (Also removing nan and Inf from label)
% (Input)   label vector (numbers, 1 x n), pattern matrix (p x n)
% (Output)  balanced new label vector and pattern matrix
% (Name, value pair) 'SamplingSize'

% [Example]
% label = [nan -1 -2 -3 Inf -1 -2 -3 Inf -1];
% exa = (1:10)';
% for i=1:10; pattern(:,i) = circshift(exa,i-1); end

% Input validity
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'label',@isvector);
addRequired(p,'pattern',@ismatrix);
addParameter(p,'SamplingSize',validScalarPosNum);
parse(p,label,pattern,varargin{:});

if iscolumn(label)
    label = label';     % label should be a row vector
end

% Validity of label & pattern matching
if length(label) ~= size(pattern,2) && length(label) ~= size(pattern,1)
    error('Mismatch(n) : label(1 x n) length & pattern(p x n) length');
elseif size(pattern,2) == size(pattern,1)
    warning('pattern is a square matrix : must be carefully used');
elseif size(pattern,1) == length(label)
    warning('label is 1 x n, but pattern is n x p : pattern is transposed');
    pattern = pattern';
end

% Removing Inf and nan
labpat1 = [label; pattern];
labpat1 = labpat1(:, all(labpat1(1:2,:)~=Inf, 1) );
labpat1 = labpat1(:, all(~isnan(labpat1(1:2,:)), 1) );
% labpat1 = labpat1(:, labpat1(1,:)~=Inf);
% labpat1 = labpat1(:, ~isnan(labpat1(1,:)));
label = labpat1(1,:);
pattern = labpat1(2:end,:);

% Determining sampling size
uniq_el = unique(label);
if ~isempty(varargin) && ismember({'SamplingSize'}, varargin{1:2:end})
    N_sampling = p.Results.SamplingSize;
else
    replab = repmat(label, length(uniq_el), 1);
    N_lab_count = sum(replab==uniq_el', 2);
    N_sampling = min(N_lab_count);  % Find smallest label number
end

% Balance the label vector by undersampling
new_labpat = nan * ones(1+size(pattern,1), N_sampling*length(uniq_el));
labpat = [label; pattern];
for el = 1:length(uniq_el)
    new_labpat(:,1+N_sampling*(el-1):N_sampling*el) = datasample(labpat(:,label==uniq_el(el)), ... 
        N_sampling, 2, 'Replace', false);
end
new_label = new_labpat(1,:);
new_pattern = new_labpat(2:end,:);

end