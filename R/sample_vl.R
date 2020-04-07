sample_vl <- function(x, ...) x[sample(length(x), ...)] # samples vectors of varying length correctly i.e. if length(x) = 1, sample won't sample from 1:x
