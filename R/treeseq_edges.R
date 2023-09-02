treeseq_edges = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_edges, ts@treeseq)
}