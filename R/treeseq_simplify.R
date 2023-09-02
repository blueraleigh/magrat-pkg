tsx_treeseq_simplify = function(ts, samples, node.map=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    storage.mode(samples) = "integer"
    handle = .Call(C_treeseq_simplify, ts@treeseq, samples)
    ts = treeseq()
    ts@treeseq = handle[[1L]]
    ts@tree = handle[[2L]]
    if (!node.map)
        return (ts)
    else
        return (structure(ts, node.map=handle[[3]]))
}
