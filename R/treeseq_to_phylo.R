treeseq_to_phylo = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    .Call(C_treeseq_to_phylo, ts@tree)
}
