#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"

/* Sample MPR histories of geographic locations of genetic ancestors
** using squared change parsimony.
**
** For each non-sample node we calculate the parameters of the quadratic 
**
**   F(x1,x2,...,xn) = p1*x1*x1 +
**                     p2*x2*x2 +
**                       ...    + 
**                     pn*xn*xn + 
**                     p1*x1    + 
**                     p2*x2    + 
**                      ...     + 
**                     pn*xn    +
**                     p
**
** Which gives the minimum sum of squared distances between all
** ancestor-descendant pairs needed to explain the spatial
** distribution of samples assuming the non-sample node is in state
** (x1, x2, ..., xn).
**
** For more information on the basic algorithm for a single tree see
**
**   Maddison, W.P. 1991. Systematic Zoology 40(3): 304-314.
**   https://www.jstor.org/stable/2992324
*/

typedef struct pdata {
    double *g;
    double *h;
    double *f;
    const double *x;
    int num_dims;
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
    int num_dims = pdata->num_dims;
    int num_dims_2 = 2*num_dims;
    int ofs = v * (num_dims_2 + 1);
    double p1, p2, b = 1;
    double *restrict h_v = pdata->h + ofs;

    if (tree->time)
    {
        b = tree->time[u] - tree->time[v];
    }
    
    if (!(tree->node_flags[v] & TSK_NODE_IS_SAMPLE))
    {
        // internal node
        double *restrict g_v = pdata->g + ofs;
    
        h_v[num_dims_2] = g_v[num_dims_2];
        for (i = 0; i < num_dims; ++i)
        {
            p1 = g_v[i];
            p2 = g_v[i+num_dims];
            h_v[i] = p1 / (b*p1 + 1);
            h_v[i+num_dims] = p2 / (b*p1 + 1);
            h_v[num_dims_2] += -(p2*p2) / (4 * (p1 + 1/b)); 
        }
    }
    else
    {
        // terminal node
        const double *restrict x = pdata->x + v*num_dims;

        h_v[num_dims_2] = 0;
        for (i = 0; i < num_dims; ++i)
        {
            h_v[i] = 1 / b;
            h_v[i+num_dims] = (-2*x[i]) / b;
            h_v[num_dims_2] += (x[i]*x[i]) / b;
        }
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    if (tree->node_flags[v] & TSK_NODE_IS_SAMPLE)
        return;

    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int num_dims = pdata->num_dims;
    int num_dims_2 = 2*num_dims;
    int num_pars = num_dims_2 + 1;
    int ofs = v * num_pars;
    tsk_id_t u = tree->parent[v];
    if (u == TSK_NULL)
    {
        memcpy(pdata->f + ofs, pdata->g + ofs, num_pars*sizeof(double));
        return;
    }
    double *restrict g_v = pdata->g + ofs;
    double *restrict h_v = pdata->h + ofs;
    double *restrict f_v = pdata->f + ofs;
    double *restrict f_u = pdata->f + u*num_pars;

    double p1;
    double p2;

    double b = 1;

    if (tree->time)
    {
        b = tree->time[u] - tree->time[v];
    }

    f_v[num_dims_2] = f_u[num_dims_2] - h_v[num_dims_2];
    for (i = 0; i < num_dims; ++i)
    {
        p1 = f_u[i] - h_v[i];
        p2 = f_u[i+num_dims] - h_v[i+num_dims];
        f_v[i] = p1 / (b*p1 + 1);
        f_v[i+num_dims] = p2 / (b*p1 + 1);
        f_v[num_dims_2] += -( p2*p2 ) / (4 * (p1 + 1/b) );
    }
    for (i = 0; i < num_pars; ++i)
        f_v[i] += g_v[i];
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int num_pars = (2*pdata->num_dims + 1);
    double *restrict g_u = pdata->g + u*num_pars;
    double *restrict h_v = pdata->h + v*num_pars;
    for (i = 0; i < num_pars; ++i)
        g_u[i] += sign*h_v[i];
}


static void
mean_tree_length(
    tsx_tree_t *tree, int t_index, double t_left, double t_right, void *params)
{
    tsk_id_t u;
    tsk_id_t virtual_root = tree->virtual_root;
    tsk_id_t *left_child = tree->left_child;
    tsk_id_t *right_sib = tree->right_sib;
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int num_dims = pdata->num_dims;
    int num_dims_2 = 2*num_dims;
    int num_pars = num_dims_2 + 1;
    double *g;
    double min_s, len = 0;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        if (!(tree->node_flags[u] & TSK_NODE_IS_SAMPLE))
        {
            g = pdata->g + u*num_pars;
            min_s = g[num_dims_2];
            for (i = 0; i < num_dims; ++i)
                min_s += -(g[i+num_dims]*g[i+num_dims]) / (4*g[i]);
            len += min_s;
        }
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

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    int i;
    int num_pars = (2*pdata->num_dims + 1);

    const double *restrict f = pdata->f + num_pars * node_id;
    double *restrict F = mpr_summary->F + num_pars * node_id;

    for (i = 0; i < num_pars; ++i)
        F[i] += w * (f[i] - F[i]) / total_weight;
}


SEXP C_treeseq_quadratic_nd_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP x)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int nprotect = 0;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);
    
    int num_dims = *INTEGER(Rf_getAttrib(x, R_DimSymbol));
    int num_pars = 2*num_dims + 1;


    SEXP g = PROTECT(Rf_allocVector(REALSXP, num_pars*num_nodes));
    SEXP h = PROTECT(Rf_allocVector(REALSXP, num_pars*num_nodes));
    SEXP f = PROTECT(Rf_allocVector(REALSXP, num_pars*num_nodes));

    SEXP F = PROTECT(Rf_allocMatrix(REALSXP, num_pars, num_nodes));
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    nprotect = 6;

    memset(REAL(g), 0, num_pars * num_nodes * sizeof(double));
    memset(REAL(h), 0, num_pars * num_nodes * sizeof(double));
    memset(REAL(f), 0, num_pars * num_nodes * sizeof(double));
    memset(REAL(F), 0, num_pars * num_nodes * sizeof(double));
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
    pdata.x = REAL(x);
    pdata.num_dims = num_dims;

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


SEXP C_treeseq_quadratic_nd_mpr_minimize(SEXP F)
{
    int i;
    int j;
    int nr = *INTEGER(Rf_getAttrib(F, R_DimSymbol));

    int num_pars = INTEGER(Rf_getAttrib(F, R_DimSymbol))[1];
    int num_dims = (num_pars - 1) / 2;

    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, nr, num_dims));
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0,
        VECTOR_ELT(Rf_getAttrib(F, R_DimNamesSymbol), 0));
    SET_VECTOR_ELT(dimnames, 1, R_NilValue);
    double *x = REAL(ans);
    const double *restrict p = REAL(F);
    for (i = 0; i < nr; ++i)
    {
        for (j = 0; j < num_dims; ++j)
            x[i+j*nr] = -p[i+(j+num_dims)*nr] / (2*p[i+j*nr]);
    }
    Rf_setAttrib(ans, R_DimNamesSymbol, dimnames);
    UNPROTECT(2);
    return ans;
}


SEXP C_treeseq_quadratic_nd_mpr_minimize_discrete(SEXP F, SEXP sites)
{
    int i;
    int j;
    int k;
    int n = *INTEGER(Rf_getAttrib(F, R_DimSymbol));
    int num_pars = INTEGER(Rf_getAttrib(F, R_DimSymbol))[1];
    int num_dims = (num_pars - 1) / 2;
    int num_sites = *INTEGER(Rf_getAttrib(sites, R_DimSymbol));
    
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, num_dims));

    double *restrict x = REAL(ans);
    const double *restrict p = REAL(F);
    const double *restrict coords = REAL(sites);

    double score;
    double min_score;

    for (i = 0; i < n; ++i)
    {
        min_score = R_PosInf;
        for (j = 0; j < num_sites; ++j)
        {
            score = p[i + (2*num_dims)*n];
            for (k = 0; k < num_dims; ++k)
            {
                score += p[i+k*n]*coords[j+k*num_sites]*coords[j+k*num_sites] +
                         p[i+(k+num_dims)*n]*coords[j+k*num_sites];
            }

            if (score < min_score)
            {
                for (k = 0; k < num_dims; ++k)
                    x[i+k*n] = coords[j+k*num_sites];
                min_score = score;
            }
        }
    }

    UNPROTECT(1);
    return ans;
}
