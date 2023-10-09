#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"


/* Sample MPR histories of geographic locations of genetic ancestors
** using generalized (Sankoff) parsimony.
**
** Geographic locations are stored as the integers 0, 1, ..., n-1, 
** where n is the number of locations.
**
** Requires a user-defined cost matrix giving the cost of transitioning
** from each state to every other.
**
** For more information on the basic algorithm for a single tree see
**
**   Sankoff, D. and P. Rousseau 1975. Math. Programming 9: 240-246.
*/

typedef struct pdata {
    int num_states;
    double *g;
    double *h;
    double *f;
    const double *cost;
} pdata_t;


typedef struct mpr_summary {
    double *F;
    double *G;
    double *node_weight;
    double *tree_length;
    double mean_tree_length;
    double tree_weight;
} mpr_summary_t;


static void
calc_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree)
{
    assert(u != TSK_NULL);
    assert(v != TSK_NULL);
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int j;
    int n = pdata->num_states;
    int vn = v*n;
    double *restrict g_v = pdata->g + vn;
    double *restrict h_v = pdata->h + vn;
    const double *restrict cost = pdata->cost;
    double t, min_t, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; i < n; ++i)
    {
        min_t = R_PosInf;
        for (j = 0; j < n; ++j)
        {
            t = scale * cost[i+j*n] + g_v[j];
            if (t < min_t)
                min_t = t;
        }
        h_v[i] = min_t;
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int j;
    tsk_id_t u = tree->parent[v];
    int n = pdata->num_states;
    int vn = v*n;
    if (u == TSK_NULL)
    {
        memcpy(pdata->f + vn, pdata->g + vn, n*sizeof(double));
        return;
    }
    int un = u*n;
    double *restrict g_v = pdata->g + vn;
    double *restrict h_v = pdata->h + vn;
    double *restrict f_v = pdata->f + vn;
    double *restrict f_u = pdata->f + un;
    const double *restrict cost = pdata->cost;
    double t, min_t, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; i < n; ++i)
    {
        min_t = R_PosInf;
        for (j = 0; j < n; ++j)
        {
            t = f_u[j] - h_v[j] + scale * cost[j+i*n] + g_v[i];
            if (t < min_t)
                min_t = t;
        }
        f_v[i] = min_t;
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int n = pdata->num_states;
    double *restrict g_u = pdata->g + u*n;
    double *restrict h_v = pdata->h + v*n;
    for (i = 0; i < n; ++i)
        g_u[i] += sign*h_v[i];
}


static void
mean_tree_length(
    tsx_tree_t *tree, int t_index, double t_left, double t_right, void *params)
{
    int i;
    tsk_id_t u;
    tsk_id_t virtual_root = tree->virtual_root;
    tsk_id_t *left_child = tree->left_child;
    tsk_id_t *right_sib = tree->right_sib;
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    double min_s;
    double *g;
    double len = 0;
    int n = pdata->num_states;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        min_s = R_PosInf;
        g = pdata->g + u*n;
        for (i = 0; i < n; ++i)
        {
            if (g[i] < min_s)
                min_s = g[i];
        }
        len += min_s;
    }
    len /= tree->num_edges;
    
    double w = t_right - t_left;

    mpr_summary_t *mpr_summary = (mpr_summary_t *)params;

    mpr_summary->tree_weight += w;
    mpr_summary->tree_length[t_index] = len;
    mpr_summary->mean_tree_length += 
        w * (len - mpr_summary->mean_tree_length) / mpr_summary->tree_weight;
}


static void
update_mean_costs(tsx_tree_t *tree, int t_index, double t_left, double t_right, 
    tsk_id_t node_id, void *params)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    mpr_summary_t *mpr_summary = (mpr_summary_t *)params;

    int i;
    int n = pdata->num_states;

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    const double *restrict f = pdata->f + n * node_id;
    const double *restrict g = pdata->g + n * node_id;
    double *restrict F = mpr_summary->F + n * node_id;
    double *restrict G = mpr_summary->G + n * node_id;

    for (i = 0; i < n; ++i)
    {
        if (R_FINITE(F[i]))
            F[i] += w * (f[i] - F[i]) / total_weight;
        if (R_FINITE(G[i]))
            G[i] += w * (g[i] - G[i]) / total_weight;
    }
}


SEXP C_treeseq_sankoff_island_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP num_states,
    SEXP g,
    SEXP cost)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int nprotect = 0;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);

    int n = *INTEGER(num_states);
    
    SEXP h = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));
    SEXP f = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));

    SEXP F = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));
    SEXP G = PROTECT(Rf_allocMatrix(REALSXP, n, num_nodes));
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    nprotect = 6;

    memset(REAL(h), 0, n * num_nodes * sizeof(double));
    memset(REAL(f), 0, n * num_nodes * sizeof(double));
    memset(REAL(F), 0, n * num_nodes * sizeof(double));
    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));
    memcpy(REAL(G), REAL(g), n * num_nodes * sizeof(double));

    mpr_summary.F = REAL(F);
    mpr_summary.G = REAL(G);
    mpr_summary.node_weight = REAL(node_weight);

    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);

    pdata.g = REAL(g);
    pdata.h = REAL(h);
    pdata.f = REAL(f);
    pdata.num_states = n;
    pdata.cost = REAL(cost);

    int rc = tsx_tree_init(&tree, ts, *INTEGER(use_brlen));

    if (rc)
        goto out;

    tree.pdata = (void *)(&pdata);
    tree.calc_stem_cost = &calc_stem_cost;
    tree.calc_final_cost = &calc_final_cost;
    tree.increment_node_cost = &increment_node_cost;

    tsx_treeseq_sankoff(
        ts,
        &tree,
        (void *)(&mpr_summary),
        &mean_tree_length,
        (void *)(&mpr_summary),
        &update_mean_costs
    );

out:
    tsx_tree_free(&tree);
    UNPROTECT(nprotect);
    return Rf_list4(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, G, F);
}


#define SQRT_DBL_EPSILON 1.490116119384765696e-8

// test for a == b
static int
fequals(double a, double b)
{
    if (!R_FINITE(a) || !R_FINITE(b)) return 0;
    // https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    double diff = fabs(a - b);
    if (diff <= SQRT_DBL_EPSILON)
    {
        return 1;
    }
    a = fabs(a);
    b = fabs(b);
    double M = (b > a) ? b : a;
    return (diff <= M * SQRT_DBL_EPSILON) ? 1 : 0;
}


// the following routines are only suitable for undirected (=symmetric) 
// state graphs

typedef struct state_graph {
    int num_states;
    // shortest path distances between all pairs of states
    double *distances;
    // adjacency matrix in compressed sparse column format
    int *i;          // row index of each non-zero weight
    int *p;          // p[j] gives the index of the weight that begins column j
    double *weights; // non-zero weights
} state_graph_t;


// sample the previous state on a shortest path from source_state to
// target_state given the current_state, which is initially set to
// the target_state
static int
prev_state(int current_state, int source_state, int target_state, 
    state_graph_t *g)
{
    int i;
    int k;
    int s;
    double d;
    double D;
    
    int num_states = g->num_states;

    // neighbors of the current state and their edge weights
    int num_neighbors = g->p[current_state+1] - g->p[current_state];
    const int *restrict neighbors = g->i + g->p[current_state];
    const double *restrict weights = g->weights + g->p[current_state];

    const double *restrict distances = g->distances;

    k = 0;

    // shortest graph distance from source state to current state
    D = distances[source_state + current_state*num_states];

    for (i = 0; i < num_neighbors; ++i)
    {
        // shortest graph distance from source state to
        // current state that goes through this neighbor
        d = distances[source_state + neighbors[i]*num_states] + weights[i];
        if (fequals(d, D))
        {
            ++k;
            if (unif_rand() < (1 / (double)k))
            {
                s = neighbors[i];
            }
        }
    }
    return s;
}


// Samples a shortest path from source to target in reverse iteration
// order. That is, to iterate from source to target use the following
// idiom:
//   
//   sample_path(source, target, ..., &path_len, path);
//
//   for (i = path_len - 1; i >= 0; --i)
//   {
//        state = path[i];
//   }
//
static void
sample_path(int source_state, int target_state, state_graph_t *g,
    int *ret_path_length, int *ret_path)
{
    int path_length = 0;
    int current_state = target_state;
    ret_path[path_length++] = current_state;
    while (current_state != source_state)
    {
        current_state = prev_state(
            current_state, source_state, target_state, g);
        ret_path[path_length++] = current_state;
    }
    *ret_path_length = path_length;
}


static void
summarize_path(int path_length, int *path, double weight, int num_states,
    double *summary_counts)
{
    int i;
    int j;
    int k;
    assert (path_length > 0);
    k = path_length - 1;
    i = path[k--];
    for (; k >= 0; --k)
    {
        j = path[k];
        summary_counts[i + j*num_states] += weight;
        i = j;
    }
}


static void
tsx_sankoff_island_sample(
    tsk_treeseq_t *ts,
    int num_samples,
    int num_nodes,
    double time_start,
    double time_end,
    const int *restrict node_states,
    state_graph_t *state_graph,
    double *ret_summary_counts)
{
    tsk_size_t i;
    int j;
    tsk_id_t u;
    tsk_id_t v;
    int to;
    int from;
    int path_length;
    int num_states = state_graph->num_states;
    int *path = malloc(num_states*sizeof(*path));
    if (!path)
    {
        TSX_WARN("memory allocation failure");
        goto out;
    }

    tsk_size_t num_edges = ts->tables->edges.num_rows;
    const tsk_id_t *restrict parent = ts->tables->edges.parent;
    const tsk_id_t *restrict child = ts->tables->edges.child;

    double branch_start;
    double branch_end;
    double branch_fraction;
    double weight;

    const double *restrict time = ts->tables->nodes.time;

    for (i = 0; i < num_edges; ++i)
    {
        u = parent[i];
        v = child[i];
        branch_start = time[v];
        branch_end = time[u];
        if (
            // +--e-------s--+
            (branch_end >= time_end && branch_start <= time_start)
            
            // +--e-------+--s
            || (branch_end >= time_end 
                && (branch_start >= time_start
                    && branch_start <= time_end))
            
            // e--+-------+--s
            || (branch_end <= time_end && branch_start >= time_start)

            // e--+-------s--+
            || ((branch_end <= time_end && branch_end >= time_start) 
                && branch_start <= time_start))
        {
            branch_fraction = (fmin2(time_end, branch_end) - 
                fmax2(time_start, branch_start)) / (branch_end - branch_start);
            weight = branch_fraction / (double)num_samples;
            for (j = 0; j < num_samples; ++j)
            {
                from = node_states[j + u*num_samples]; 
                to = node_states[j + v*num_samples];
                sample_path(from, to, state_graph, &path_length, path);
                summarize_path(path_length, path, weight, num_states, 
                    ret_summary_counts);
            }
        }
    }
out:
    free(path);
}


SEXP C_treeseq_sankoff_island_mpr_sample(
    SEXP treeseq,
    SEXP node_states,
    SEXP time_start,
    SEXP time_end,
    SEXP cost,
    SEXP A)
{
    state_graph_t graph;
    
    int num_samples = INTEGER(Rf_getAttrib(node_states, R_DimSymbol))[0];
    int num_nodes = INTEGER(Rf_getAttrib(node_states, R_DimSymbol))[1];
    int num_states = *INTEGER(Rf_getAttrib(cost, R_DimSymbol));

    graph.num_states = num_states;
    graph.distances = REAL(cost);
    graph.i = INTEGER(R_do_slot(A, Rf_mkString("i")));
    graph.p = INTEGER(R_do_slot(A, Rf_mkString("p")));
    graph.weights = REAL(R_do_slot(A, Rf_mkString("x")));

    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);

    SEXP ret = PROTECT(Rf_allocMatrix(REALSXP, num_states, num_states));

    double *summary_counts = REAL(ret);
    memset(summary_counts, 0, num_states*num_states*sizeof(double));

    GetRNGstate();
    tsx_sankoff_island_sample(
        ts,
        num_samples,
        num_nodes,
        *REAL(time_start),
        *REAL(time_end),
        INTEGER(node_states),
        &graph,
        summary_counts
    );
    PutRNGstate();
    
    UNPROTECT(1);
    return ret;
}
