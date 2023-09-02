treeseq_sample = function(ts, at=-1L)
{
    stopifnot(inherits(ts, "treeseq"))
    storage.mode(at) = "integer"
    
    .Call(C_treeseq_sample, ts@tree, at)

    invisible()
}
