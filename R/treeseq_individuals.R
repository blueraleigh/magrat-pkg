treeseq_individuals = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_individuals, ts@treeseq)
}