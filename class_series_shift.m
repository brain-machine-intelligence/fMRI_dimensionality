function shifted_series = class_series_shift(class_series, varargin)

chng_idx = [find([1 diff(class_series)~=0]~=0) length(class_series)]; % length = # of block change + 1
blkLens = diff(chng_idx);   % each block length
shifted_series = nan(size(class_series));
for i = 1:(length(blkLens) - 1)
    curr_class = class_series(chng_idx(i));
    next_changed = chng_idx(i+1);
    next_blkLen = blkLens(i+1);
    shifted_series(next_changed:(next_changed+next_blkLen-1)) = curr_class;
end


