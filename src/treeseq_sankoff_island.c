#include <R.h>
#include <Rinternals.h>
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
    double *restrict F = mpr_summary->F + n * node_id;

    for (i = 0; i < n; ++i)
    {
        if (R_FINITE(F[i]))
            F[i] += w * (f[i] - F[i]) / total_weight;
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
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    nprotect = 5;

    memset(REAL(h), 0, n * num_nodes * sizeof(double));
    memset(REAL(f), 0, n * num_nodes * sizeof(double));
    memset(REAL(F), 0, n * num_nodes * sizeof(double));
    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.F = REAL(F);
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
    return Rf_list3(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, F);
}
