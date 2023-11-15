interleave_matrices <- function(L){
    return(do.call(rbind, L)[order(sequence(sapply(L, nrow))),])
}
