treeseq_nodes = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_nodes, ts@treeseq)
}