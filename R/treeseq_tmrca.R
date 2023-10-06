treeseq_tmrca = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_tmrca, ts@tree)
}