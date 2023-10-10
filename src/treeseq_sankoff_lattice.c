#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"


/* Sample MPR histories of geographic locations of genetic ancestors
** using generalized (Sankoff) parsimony.
**
** Geographic locations are stored as row and column indices of
** a rectangular grid representing discretized geographic space.
**
** In a single dimension, dispersal events occur between
** neighboring states according to the following diagram:
**
** [ 0 ] <-> [ 1 ] <-> [ 2 ] <-> ... <-> [ n-1 ]
**
** States 0 and n-1 communicate when boundaries are periodic.
**
** Cost matrices are not stored explictly but can be worked out from
** the indices of the geographic states as we assign a cost of 1 for
** each dispersal.
**
** For more information on the basic algorithm for a single tree see
**
**   Sankoff, D. and P. Rousseau 1975. Math. Programming 9: 240-246.
*/


typedef struct pdata {
    int num_states_x;
    int num_states_y;
    int periodic_x;
    int periodic_y;
    double *gx;
    double *gy;
    double *hx;
    double *hy;
    double *fx;
    double *fy;
} pdata_t;


typedef struct mpr_summary {
    double *Fx;
    double *Fy;
    double *node_weight;
    double *tree_length;
    double mean_tree_length;
    double tree_weight;
} mpr_summary_t;


static double
cost(int i, int j, int n, int periodic)
{
    int d = abs(j - i);
    if (periodic)
    {
        int h = n / 2;
        return (d <= h) ? (double)d : (double)(n - d);
    }
    return (double)d;
}


static void
calc_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree)
{
    assert(u != TSK_NULL);
    assert(v != TSK_NULL);
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int j;
    int nx = pdata->num_states_x;
    int ny = pdata->num_states_y;
    int px = pdata->periodic_x;
    int py = pdata->periodic_y;
    int vnx = v*nx;
    int vny = v*ny;
    double *restrict gx_v = pdata->gx + vnx;
    double *restrict gy_v = pdata->gy + vny;
    double *restrict hx_v = pdata->hx + vnx;
    double *restrict hy_v = pdata->hy + vny;
    double tx, ty, min_tx, min_ty, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; (i < nx) && (i < ny); ++i)
    {
        min_tx = R_PosInf;
        min_ty = R_PosInf;
        for (j = 0; (j < nx) && (j < ny); ++j)
        {
            tx = scale * cost(i, j, nx, px) + gx_v[j];
            ty = scale * cost(i, j, ny, py) + gy_v[j];
            if (tx < min_tx)
                min_tx = tx;
            if (ty < min_ty)
                min_ty = ty;
        }
        for (; j < nx; ++j)
        {
            tx = scale * cost(i, j, nx, px) + gx_v[j];
            if (tx < min_tx)
                min_tx = tx;
        }
        for (; j < ny; ++j)
        {
            ty = scale * cost(i, j, ny, py) + gy_v[j];
            if (ty < min_ty)
                min_ty = ty;
        }
        hx_v[i] = min_tx;
        hy_v[i] = min_ty;
    }
    for (; i < nx; ++i)
    {
        min_tx = R_PosInf;
        for (j = 0; j < nx; ++j)
        {
            tx = scale * cost(i, j, nx, px) + gx_v[j];
            if (tx < min_tx)
                min_tx = tx;
        }
        hx_v[i] = min_tx;
    }
    for (; i < ny; ++i)
    {
        min_ty = R_PosInf;
        for (j = 0; j < ny; ++j)
        {
            ty = scale * cost(i, j, ny, py) + gy_v[j];
            if (ty < min_ty)
                min_ty = ty;
        }
        hy_v[i] = min_ty;
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    int i;
    int j;
    tsk_id_t u = tree->parent[v];
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int nx = pdata->num_states_x;
    int ny = pdata->num_states_y;
    int px = pdata->periodic_x;
    int py = pdata->periodic_y;
    int vnx = v*nx;
    int vny = v*ny;
    if (u == TSK_NULL)
    {
        memcpy(pdata->fx + vnx, pdata->gx + vnx, nx*sizeof(double));
        memcpy(pdata->fy + vny, pdata->gy + vny, ny*sizeof(double));
        return;
    }
    int unx = u*nx;
    int uny = u*ny;
    double *restrict gx_v = pdata->gx + vnx;
    double *restrict gy_v = pdata->gy + vny;
    double *restrict hx_v = pdata->hx + vnx;
    double *restrict hy_v = pdata->hy + vny;
    double *restrict fx_v = pdata->fx + vnx;
    double *restrict fy_v = pdata->fy + vny;
    double *restrict fx_u = pdata->fx + unx;
    double *restrict fy_u = pdata->fy + uny;
    double tx, ty, min_tx, min_ty, scale = 1;
    if (tree->time)
    {
        scale = 1 / (tree->time[u] - tree->time[v]);
    }
    for (i = 0; (i < nx) && (i < ny); ++i)
    {
        min_tx = R_PosInf;
        min_ty = R_PosInf;
        for (j = 0; (j < nx) && (j < ny); ++j)
        {
            tx = fx_u[j] - hx_v[j] + scale * cost(j, i, nx, px) + gx_v[i];
            ty = fy_u[j] - hy_v[j] + scale * cost(j, i, ny, py) + gy_v[i];
            if (tx < min_tx)
                min_tx = tx;
            if (ty < min_ty)
                min_ty = ty;
        }
        for (; j < nx; ++j)
        {
            tx = fx_u[j] - hx_v[j] + scale * cost(j, i, nx, px) + gx_v[i];
            if (tx < min_tx)
                min_tx = tx;
        }
        for (; j < ny; ++j)
        {
            ty = fy_u[j] - hy_v[j] + scale * cost(j, i, ny, py) + gy_v[i];
            if (ty < min_ty)
                min_ty = ty;
        }
        fx_v[i] = min_tx;
        fy_v[i] = min_ty;
    }
    for (; i < nx; ++i)
    {
        min_tx = R_PosInf;
        for (j = 0; j < nx; ++j)
        {
            tx = fx_u[j] - hx_v[j] + scale * cost(j, i, nx, px) + gx_v[i];
            if (tx < min_tx)
                min_tx = tx;
        }
        fx_v[i] = min_tx;
    }
    for (; i < ny; ++i)
    {
        min_ty = R_PosInf;
        for (j = 0; j < ny; ++j)
        {
            ty = fy_u[j] - hy_v[j] + scale * cost(j, i, ny, py) + gy_v[i];
            if (ty < min_ty)
                min_ty = ty;
        }
        fy_v[i] = min_ty;
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int i;
    int nx = pdata->num_states_x;
    int ny = pdata->num_states_y;
    int unx = u*nx;
    int uny = u*ny;
    int vnx = v*nx;
    int vny = v*ny;
    double *restrict gx_u = pdata->gx + unx;
    double *restrict gy_u = pdata->gy + uny;
    double *restrict hx_v = pdata->hx + vnx;
    double *restrict hy_v = pdata->hy + vny;
    for (i = 0; (i < nx) && (i < ny); ++i)
    {
        gx_u[i] += sign*hx_v[i];
        gy_u[i] += sign*hy_v[i];
    }
    for (; i < nx; ++i)
        gx_u[i] += sign*hx_v[i];
    for (; i < ny; ++i)
        gy_u[i] += sign*hy_v[i];
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
    double minx;
    double miny;
    double *gx;
    double *gy;
    double len = 0;
    int nx = pdata->num_states_x;
    int ny = pdata->num_states_y;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        minx = R_PosInf;
        miny = R_PosInf;
        gx = pdata->gx + u*nx;
        gy = pdata->gy + u*ny;
        for (i = 0; i < nx && i < ny; ++i)
        {
            if (gx[i] < minx)
                minx = gx[i];
            if (gy[i] < miny)
                miny = gy[i];
        }
        for (; i < nx; ++i)
        {
            if (gx[i] < minx)
                minx = gx[i];
        }
        for (; i < ny; ++i)
        {
            if (gy[i] < miny)
                miny = gy[i];
        }
        len += minx + miny;
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
    int nx = pdata->num_states_x;
    int ny = pdata->num_states_y;

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    const double *restrict fx = pdata->fx + nx * node_id;
    const double *restrict fy = pdata->fy + ny * node_id;
    double *restrict Fx = mpr_summary->Fx + nx * node_id;
    double *restrict Fy = mpr_summary->Fy + ny * node_id;

    for (i = 0; i < nx && i < ny; ++i)
    {
        if (R_FINITE(Fx[i]))
            Fx[i] += w * (fx[i] - Fx[i]) / total_weight;
        if (R_FINITE(Fy[i]))
            Fy[i] += w * (fy[i] - Fy[i]) / total_weight;
    }
    for (; i < nx; ++i)
    {
        if (R_FINITE(Fx[i]))
            Fx[i] += w * (fx[i] - Fx[i]) / total_weight;
    }
    for (; i < ny; ++i)
    {
        if (R_FINITE(Fy[i]))
            Fy[i] += w * (fy[i] - Fy[i]) / total_weight;
    }
}


SEXP C_treeseq_lattice_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP num_states_x,
    SEXP num_states_y,
    SEXP px,
    SEXP py,
    SEXP gx,
    SEXP gy)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int nprotect = 0;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);

    int nx = *INTEGER(num_states_x);
    int ny = *INTEGER(num_states_y);
    
    SEXP hx = PROTECT(Rf_allocMatrix(REALSXP, nx, num_nodes));
    SEXP hy = PROTECT(Rf_allocMatrix(REALSXP, ny, num_nodes));
    SEXP fx = PROTECT(Rf_allocMatrix(REALSXP, nx, num_nodes));
    SEXP fy = PROTECT(Rf_allocMatrix(REALSXP, ny, num_nodes));

    SEXP Fx = PROTECT(Rf_allocMatrix(REALSXP, nx, num_nodes));
    SEXP Fy = PROTECT(Rf_allocMatrix(REALSXP, ny, num_nodes));
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    nprotect = 8;

    memset(REAL(hx), 0, nx * num_nodes * sizeof(double));
    memset(REAL(hy), 0, ny * num_nodes * sizeof(double));
    memset(REAL(fx), 0, nx * num_nodes * sizeof(double));
    memset(REAL(fy), 0, ny * num_nodes * sizeof(double));
    memset(REAL(Fx), 0, nx * num_nodes * sizeof(double));
    memset(REAL(Fy), 0, ny * num_nodes * sizeof(double));
    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.Fx = REAL(Fx);
    mpr_summary.Fy = REAL(Fy);
    mpr_summary.node_weight = REAL(node_weight);

    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);

    pdata.gx = REAL(gx);
    pdata.gy = REAL(gy);
    pdata.hx = REAL(hx);
    pdata.hy = REAL(hy);
    pdata.fx = REAL(fx);
    pdata.fy = REAL(fy);
    pdata.num_states_x = nx;
    pdata.num_states_y = ny;
    pdata.periodic_x = *INTEGER(px);
    pdata.periodic_y = *INTEGER(py);

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
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, Fx, Fy);
}
