function [new_label, newlab_info] = bin_labels_combination(bin_labels_mat, varargin)
%Given m binary label vectors (values : 0, 1), make a single new label vector based on their combinatory pattern (values : 1, 2, ... , 2^m)
% (input) bin_labels_mat : m X N_timepoints(N_slices) (values = 0, 1)
% (output) new_label : 1 X N_timepoints (possible values = 1, 2, ... , 2^m, null label = NaN)
    
    % Input validity
    p = inputParser;
    validScalarNum = @(x) isnumeric(x) && isscalar(x);
    addRequired(p,'bin_labels_mat');
    addParameter(p,'Verbose',validScalarNum);
    parse(p,bin_labels_mat,varargin{:});
    
    % Setting default values for options
    options = struct('Verbose',1);
    optionNames = fieldnames(options);
    
    % varargin validity
    for pair = reshape(varargin, 2, [])
        if ismember(pair{1}, optionNames)
            options.(pair{1}) = pair{2};
        else
            error('(bin_labels_combination) %s is not recognized parameter name', pair{1});
        end
    end

    if ~all(ismember(unique(bin_labels_mat(~isnan(bin_labels_mat))),[0,1]))
        disp(['(bin_labels_combination) Warning: bin_labels_mat has ', num2str(unique(bin_labels_mat))]);
    end
    
    m = size(bin_labels_mat, 1);
    N_timepoints = size(bin_labels_mat, 2);
    bin_con = cell(1,m);
    for i=1:m
        bin_con{i} = [1 0];
    end
    label_patterns = combvec(bin_con{:});       % m X 2^m
    
    new_label = nan * ones(1,N_timepoints);
    newlab_info.label_combinations = zeros(size(label_patterns));
    newlab_info.newlab_comb_map = zeros(size(label_patterns,1)+1, size(label_patterns,2));
    for ii=1:size(label_patterns,2)
        ii_labelled_idx = logical(prod(bin_labels_mat==label_patterns(:,ii), 1));
        new_label(ii_labelled_idx) = ii;
        newlab_info.label_combinations(:,ii) = label_patterns(:,ii);
        newlab_info.newlab_comb_map(1,ii) = ii;
        newlab_info.newlab_comb_map(2:end,ii) = label_patterns(:,ii);
    end
    
    if options.('Verbose')
        disp('(bin_labels_combination) new_label & bin_labels mapping (first row : new labels, below matrix : prior label combinations) ');
        disp(num2str(newlab_info.newlab_comb_map));
    end
end

