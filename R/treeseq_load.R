treeseq_load = function(filename)
{
    handle = .Call(C_treeseq_load, filename)
    ts = treeseq()
    ts@treeseq = handle[[1L]]
    ts@tree = handle[[2L]]
    ts
}
