treeseq_intervals = function(ts)
{
    stopifnot(inherits(ts, "treeseq"))
    structure(
        .Call(C_treeseq_intervals, ts@tree)
        , dimnames=list(
            NULL
            , c("left","right","length","num_edges","num_roots","max_root_age")
        )
    )
}
