treeseq_sankoff_island_mpr = function(ts, sample_locations, cost, 
    use_brlen=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(is.matrix(cost))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(all(colnames(sample_locations) %in% c("node_id","state_id")))
    storage.mode(sample_locations) = "integer"
    stopifnot(nrow(cost) == ncol(cost))
    stopifnot(all(cost >= 0L))
    num_states = nrow(cost)
    stopifnot(all(sample_locations[,"state_id"] > 0L))
    stopifnot(all(sample_locations[,"state_id"] <= num_states))
    N = nrow(treeseq_nodes(ts))
    G = matrix(0, num_states, N)
    sample_ids = sample_locations[, "node_id"] + 1L
    for (i in 1:nrow(sample_locations))
    {
        k = sample_ids[i]
        j = sample_locations[i, "state_id"]
        G[, k] = Inf
        G[j, k] = 0
    }
    L = .Call(
        C_treeseq_sankoff_island_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , num_states,
        , G
        , cost
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr")
    L[[3]] = t(L[[3]][, -sample_ids])
    rownames(L[[3]]) = (0:(N-1))[-sample_ids]
    L
}


treeseq_sankoff_planar_mpr = function(ts, sample_locations,
    use_brlen=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(all(colnames(sample_locations) %in% c("node_id","x", "y")))
    storage.mode(sample_locations) = "double"
    N = nrow(treeseq_nodes(ts))
    num_samples = nrow(sample_locations)
    x = numeric(N)
    y = numeric(N)
    sample_ids = sample_locations[, "node_id"] + 1L
    x[sample_ids] = sample_locations[,"x"]
    y[sample_ids] = sample_locations[,"y"]
    L = .Call(
        C_treeseq_sankoff_planar_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , x
        , y
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr")
    L[[3]] = t(L[[3]][, -sample_ids])
    rownames(L[[3]]) = (0:(N-1))[-sample_ids]
    L
}


treeseq_sankoff_lattice_mpr = function(
    ts,
    num_states_x,
    num_states_y,
    periodic_x,
    periodic_y,
    sample_locations,
    use_brlen=FALSE)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(sample_locations))
    stopifnot(!is.null(colnames(sample_locations)))
    stopifnot(all(colnames(sample_locations) %in%
        c("node_id","xmin","xmax","ymin","ymax")))
    storage.mode(sample_locations) = "integer"
    stopifnot(all(sample_locations[,"xmin"] > 0L & 
        sample_locations[,"xmin"] <= num_states_x))
    stopifnot(all(sample_locations[,"xmax"] > 0L & 
        sample_locations[,"xmax"] <= num_states_x))
    stopifnot(all(sample_locations[,"ymin"] > 0L & 
        sample_locations[,"ymin"] <= num_states_y))
    stopifnot(all(sample_locations[,"ymax"] > 0L & 
        sample_locations[,"ymax"] <= num_states_y))
    stopifnot(num_states_x > 1)
    stopifnot(num_states_y > 1)
    N = nrow(treeseq_nodes(ts))
    num_samples = nrow(sample_locations)
    Gx = matrix(0, num_states_x, N)
    Gy = matrix(0, num_states_y, N)
    sample_ids = sample_locations[, "node_id"] + 1L
    for (i in 1:num_samples)
    {
        k = sample_ids[i]
        j = sample_locations[i,"xmin"]:sample_locations[i,"xmax"]
        Gx[, k] = Inf
        Gx[j, k] = 0
        j = sample_locations[i,"ymin"]:sample_locations[i,"ymax"]
        Gy[, k] = Inf
        Gy[j, k] = 0
    }
    storage.mode(num_states_x) = "integer"
    storage.mode(num_states_y) = "integer"
    storage.mode(periodic_x) = "integer"
    storage.mode(periodic_y) = "integer"
    L = .Call(
        C_treeseq_sankoff_lattice_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , num_states_x
        , num_states_y
        , periodic_x
        , periodic_y
        , Gx
        , Gy
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr_x", "mpr_y")
    L[[3]] = t(L[[3]][, -sample_ids])
    L[[4]] = t(L[[4]][, -sample_ids])
    rownames(L[[3]]) = rownames(L[[4]]) = (0:(N-1))[-sample_ids]
    L
}
