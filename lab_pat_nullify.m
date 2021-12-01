function pattern = lab_pat_nullify(label, pattern, varargin)
% nullify input label's effect in the given patterns
% label vector and pattern matrix. (Also removing nan and Inf from label)
% (Input)   discrete label vector (numbers, 1 x n), pattern matrix or tensor
% (Output)  nullified pattern matrix
% (Name, value pair) 'NullDim'

% [Example]
% label = [nan -1 -2 -3 Inf -1 -2 -3 Inf -1];
% exa = (1:10)';
% for i=1:10; pattern(:,i) = circshift(exa,i-1); end

% Input validity
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p,'label',@isvector);
% addRequired(p,'pattern',@ismatrix);
addParameter(p,'NullDim',validScalarPosNum);
parse(p,label,varargin{:});

if ~isempty(varargin) && ismember({'NullDim'}, varargin{1:2:end})
    nulldim = p.Results.NullDim;
else
    nulldim = 2;
end

if iscolumn(label)
    label = label';     % label should be a row vector
end

% Validity of label & pattern matching
if ~any(ismember(size(pattern), length(label)))
    error('Mismatch(n) : label(1 x n) length & pattern(p x n) length');
elseif size(pattern,2) == size(pattern,1)
    warning('pattern is a square matrix : must be carefully used');
elseif size(pattern,1) == length(label)
    warning('label is 1 x n, but pattern is n x p : pattern is transposed');
    pattern = pattern';
end

% Removing Inf and nan
switch nulldim
    case 2
%         pattern(:, any(isnan(pattern(1:2, :)), 1)) = [];
    case 3
%         error('(lab_pat_nullify) not implemented')
end

% nullifying (mean-centering per label class)
[unqSet, ~] = unq_elms(label);
for ui = 1:length(unqSet)
    if nulldim == 2
        target_pat = pattern(:, label == unqSet(ui));
%         pattern(:, label == unqSet(ui)) = target_pat - mean(target_pat, nulldim);
        pattern(:, label == unqSet(ui)) = target_pat - nanmean(target_pat, nulldim);
    elseif nulldim == 3
        target_pat = pattern(:, :, label == unqSet(ui));
%         pattern(:, :, label == unqSet(ui)) = target_pat - mean(target_pat, nulldim);
        pattern(:, :, label == unqSet(ui)) = target_pat - nanmean(target_pat, nulldim);
    else
        error('(lab_pat_nullify) NullDim must be 2 or 3')
    end
end

