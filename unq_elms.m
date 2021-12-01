function [unqSet, unqNs] = unq_elms(targ_vec, varargin)
% Finding and counting unique elements in the target vector

% Input validity
p = inputParser;
validScalarNum = @(x) isnumeric(x) && isscalar(x);
addRequired(p,'targ_vec',@isvector);
addParameter(p,'NaNout',validScalarNum);
addParameter(p,'Verbose',validScalarNum);
parse(p,targ_vec,varargin{:});

% Setting default values for options
options = struct('NaNout',1,'Verbose',0);
optionNames = fieldnames(options);

n_args = length(varargin);
if round(n_args/2) ~= n_args/2
    fprintf('n_args=%d \n',n_args)
    error('(unq_elms) varargin pair error (NaNout pair & Verbose pair)')
end

% 
% if ~isempty(varargin) && ismember({'NaNout'}, varargin{1:2:end})
%     nanout = p.Results.NaNout;
% else
%     nanout = 1;
% end

for pair = reshape(varargin, 2, [])
    if ismember(pair{1}, optionNames)
        options.(pair{1}) = pair{2};
    else
        error('(unq_elms) %s is not recognized parameter name', pair{1});
    end
end
    
% 
% if ~isempty(varargin) && ismember({'Verbose'}, varargin{1:2:end})
%     printout = p.Results.Verbose;
% else
%     printout = 0;
% end

if options.('NaNout')
    [unqSet,~,idx] = unique(targ_vec(~isnan(targ_vec)));
else
    [unqSet,~,idx] = unique(targ_vec);
end

unqNs = accumarray(idx,1);

if options.('Verbose')
    for el = 1:length(unqSet)
        if iscell(unqSet)
            fprintf([num2str(unqSet{el}) ' ' num2str(unqNs(el)) '\n']);
        else
            fprintf([num2str(unqSet(el)) ' ' num2str(unqNs(el)) '\n']);
        end
    end
end

