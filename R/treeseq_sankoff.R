treeseq_island_mpr = function(ts, sample_locations, cost, 
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
    storage.mode(cost) = "double"
    num_states = nrow(cost)
    stopifnot(all(sample_locations[,"state_id"] > 0L))
    stopifnot(all(sample_locations[,"state_id"] <= num_states))
    N = nrow(treeseq_nodes(ts))
    G = matrix(0, num_states, N)
    sample_ids = sample_locations[, "node_id"] + 1L
    state_ids = sample_locations[, "state_id"]
    G[, sample_ids] = Inf
    G[cbind(state_ids, sample_ids)] = 0
    structure(.Call(
        C_treeseq_island_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , num_states
        , G
        , cost
    ), names=c("mean_tree_length", "tree_length", "G", "F"))
}


treeseq_island_mpr_sample = function(ts, F, cost, adjacency_matrix,
    num_samples=1L, time_start=0, time_end=Inf)
{
    stopifnot(inherits(ts, "treeseq"))
    stopifnot(is.matrix(cost))
    stopifnot(is.matrix(adjacency_matrix))
    stopifnot(nrow(cost) == ncol(cost))
    stopifnot(all(cost >= 0L))
    stopifnot(all(rowSums(adjacency_matrix) > 0))
    stopifnot(all(colSums(adjacency_matrix) > 0))
    stopifnot(time_start < time_end)
    stopifnot(time_start >= 0)
    stopifnot(time_end >= 0)
    storage.mode(cost) = "double"
    storage.mode(adjacency_matrix) = "double"
    A = Matrix::Matrix(adjacency_matrix, sparse=TRUE)
    stopifnot(inherits(A, "dgCMatrix"))

    num_states = nrow(adjacency_matrix)
    state_space = 0:(num_states - 1)
    state_weights = apply(F, 2L, function(s) {
        p = exp(-s + min(s))
        p / sum(p)
    })

    node_states = apply(state_weights, 2L, function(p) {
        sample(state_space, num_samples, replace=TRUE, prob=p)
    })

    if (num_samples == 1L)
        node_states = t(node_states)

    .Call(
        C_treeseq_island_mpr_sample
        , ts@treeseq
        , node_states
        , time_start
        , time_end
        , cost
        , A
    )
}


treeseq_lattice_mpr = function(
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
        C_treeseq_lattice_mpr
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


treeseq_quadratic_mpr = function(ts, sample_locations, use_brlen=FALSE)
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
        C_treeseq_quadratic_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , x
        , y
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr")
    L[[3]] = t(L[[3]][, -sample_ids])
    rownames(L[[3]]) = (0:(N-1))[-sample_ids]
    structure(L, class=c("quadratic", "mpr"))
}


treeseq_quadratic_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "quadratic") && inherits(obj, "mpr"))
    .Call(C_treeseq_quadratic_mpr_minimize, obj$mpr)
}

treeseq_quadratic_mpr_sample = function(obj, rate)
{
    stopifnot(inherits(obj, "quadratic") && inherits(obj, "mpr"))
    if (missing(rate))
        rate = 1
    .Call(C_treeseq_quadratic_mpr_sample, obj$mpr, rate)
}


treeseq_linear_mpr = function(ts, sample_locations, use_brlen=FALSE)
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
    nx = length(unique(sample_locations[,"x"]))
    ny = length(unique(sample_locations[,"y"]))
    L = .Call(
        C_treeseq_linear_mpr
        , ts@treeseq
        , as.integer(use_brlen)
        , x
        , y
        , nx
        , ny
    )
    names(L) = c("mean_tree_length", "tree_length", "mpr_x", "mpr_y")
    structure(L, class=c("linear", "mpr"))
}


treeseq_linear_mpr_minimize = function(obj)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    .Call(C_treeseq_linear_mpr_minimize, obj$mpr_x, obj$mpr_y)
}


treeseq_linear_mpr_sample = function(obj, rate)
{
    stopifnot(inherits(obj, "linear") && inherits(obj, "mpr"))
    if (missing(rate))
        rate = 1
    .Call(C_treeseq_linear_mpr_sample, obj$mpr_x, obj$mpr_y, rate)
}
