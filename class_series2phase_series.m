function phase_series = class_series2phase_series(class_series, varargin)

options = struct('PhaseNumber', 2);
option_names = fieldnames(options);
if mode(length(varargin), 2) == 1
    error('(class_series2phase_series) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(class_series2phase_series) %s is not a recognized parameter name', pair{1})
    end
end
n_phase = options.PhaseNumber;

chng_idx = [find([1 diff(class_series)~=0]~=0) length(class_series)]; % length = # of block change + 1
blkLens = diff(chng_idx);   % each block length
phase_series = n_phase * ones(size(class_series)); % 0:early, 1:late
%             blkEarly = 2 * ones(size(blk_con)); % 1:early(33%), 2:late(66%)
%             blk3phase = 3 * ones(size(blk_con));
for i = 1:length(blkLens)
    changed = chng_idx(i);
    blkLen = blkLens(i);
    for j = 1:n_phase
        start = changed + (j - 1) * ceil(blkLen / n_phase);
        end_ = min(changed + blkLen - 1, start + ceil(blkLen / n_phase) - 1);
        phase_series(start:end_) = j;
    end
end


