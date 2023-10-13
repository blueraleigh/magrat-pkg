treeseq_populations = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_populations, ts@treeseq)
}
