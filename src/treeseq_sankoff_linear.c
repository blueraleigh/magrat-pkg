#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <tskit.h>

#include "treeseq_sankoff.h"
#include "error.h"
#include "fequals.h"

/* Sample MPR histories of geographic locations of genetic ancestors
** using linear parsimony.
**
** For each non-sample node we calculate the parameters of a convex
** piecewise linear function 


       { s[0] + a[0] * x                 if x <= x[1]
       { s[1] + a[1] * (x - x[1])        if x[1] < x <= x[2] 
F(x) = { ...
       { s[k-1] + a[k-1] * (x - x[k-1])  if x[k-1] < x <= x[k]
       { s[k]   + a[k]   * (x - x[k])    if x > x[k]

where s[i+1] = s[i] + a[i]*(x[i+1] - x[i])

** Which gives the minimum sum of distances between all
** ancestor-descendant pairs needed to explain the spatial
** distribution of samples assuming the non-sample node is in
** location x.
**
** For more information on the basic algorithm for a single tree see
**
**  Miklós Csurös. 2008. Ancestral Reconstruction by Asymmetric Wagner
**  Parsimony over Continuous Characters and Squared Parsimony over
**  Distributions. Pp 72-86. In: Nelson, C.E., Vialette, S. (eds)
**  Comparative Genomics.
*/


// convex piecewise linear function.
// slopes and breakpoints are stored in non-decreasing order
typedef struct plf {
    double *slope;
    double *breakpoint;
    double intercept;
    int num_breaks;
} plf_t;


typedef struct pdata {
    plf_t *gx;
    plf_t *hx;
    plf_t *fx;
    plf_t *gy;
    plf_t *hy;
    plf_t *fy;
    // workspace
    plf_t *tmp;
    // observed states
    const double *x;
    const double *y;
} pdata_t;


typedef struct mpr_summary {
    plf_t *Fx;
    plf_t *Fy;
    plf_t *tmp;
    double *node_weight;
    double *tree_length;
    double mean_tree_length;
    double tree_weight;
} mpr_summary_t;


static int
plf_init(plf_t *f, int n)
{
    f->slope = calloc(n+1, sizeof(double));
    f->breakpoint = calloc(n, sizeof(double));
    f->intercept = 0;
    f->num_breaks = 0;
    if (f->slope == NULL || f->breakpoint == NULL)
    {
        TSX_WARN("memory allocation failed");
        return 1;
    }
    return 0;
}


static void
plf_free(plf_t *f)
{
    free(f->slope);
    free(f->breakpoint);
}


static double
plf_min_score(plf_t *f)
{
    int i;
    double shift;
    double *restrict fs = f->slope;
    double *restrict fb = f->breakpoint;
    i = 0;
    shift = f->intercept + fs[0]*fb[0];
    while (fs[i+1] < 0)
    {
        ++i;
        shift += fs[i] * (fb[i] - fb[i-1]);
    }
    return shift;
}


// ret(x) = min_z { c(x,z) + f(z) } where the cost c(x,z)
// of going from x to z is b * abs(x - z)
static void
plf_min(plf_t *f, double b, plf_t *ret)
{
    assert(b > 0);
    assert(f != ret);
    int i;
    int x_left;
    int num_breaks = 0;
    double shift;
    double minus_b = -b;

    double *restrict fs = f->slope;
    double *restrict fb = f->breakpoint;

    double *restrict rs = ret->slope;
    double *restrict rb = ret->breakpoint;

    if (minus_b <= fs[0])
    {
        x_left = 0;
        ret->intercept = f->intercept;
    }
    else
    {
        i = 0;
        shift = f->intercept + fs[0]*fb[0];
        while (minus_b >= fs[i+1])
        {
            ++i;
            shift += fs[i] * (fb[i] - fb[i-1]);
        }
        x_left = i + 1;
        ret->intercept = shift + b * fb[i];
        rs[0] = minus_b;
        rb[0] = fb[i];
        ++num_breaks;
    }
    if (b >= fs[f->num_breaks])
    {
        for (i = x_left; i < f->num_breaks; ++i)
        {
            rs[num_breaks] = fs[i];
            rb[num_breaks] = fb[i];
            ++num_breaks;
        }
        rs[num_breaks] = fs[f->num_breaks];
    }
    else
    {
        i = x_left;
        while (b >= fs[i])
        {
            assert(i < f->num_breaks);
            rs[num_breaks] = fs[i];
            rb[num_breaks] = fb[i];
            ++num_breaks;
            ++i;
        }
        // check for colinearity to avoid superfluous breakpoints
        if (fequals(b, rs[num_breaks - 1]))
            --num_breaks;
        else
            rs[num_breaks] = b;
    }
    ret->num_breaks = num_breaks;
}


// ret = a + sign*scale*b
static void
plf_add(plf_t *a, plf_t *b, int sign, double scale, plf_t *ret)
{
    assert(a != ret);
    assert(b != ret);
    assert(a != b);
    int i;
    int j;
    int num_breaks;

    double prev_slope;
    double next_slope;
    double next_break;

    double *restrict breakpoint = ret->breakpoint;
    double *restrict slope = ret->slope;

    double *restrict ab = a->breakpoint;
    double *restrict as = a->slope;
    double *restrict bb = b->breakpoint;
    double *restrict bs = b->slope;

    i = 0;
    j = 0;
    num_breaks = 0;
    prev_slope = R_NegInf;

    while (i < a->num_breaks && j < b->num_breaks)
    {
        if (ab[i] < bb[j])
        {
            next_break = ab[i];
            next_slope = as[i] + sign * scale * bs[j];
            ++i;
        }
        else if (ab[i] > bb[j])
        {
            next_break = bb[j];
            next_slope = as[i] + sign * scale * bs[j];
            ++j;
        }
        else
        {
            next_break = ab[i];
            next_slope = as[i] + sign * scale * bs[j];
            ++i;
            ++j;
        }
        if (fequals(prev_slope, next_slope))
        {
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
        }
    }
    
    while (i < a->num_breaks)
    {
        next_break = ab[i];
        next_slope = as[i] + sign * scale * bs[j];
        ++i;
        if (fequals(prev_slope, next_slope))
        {
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
        }
    }
    while (j < b->num_breaks)
    {
        next_break = bb[j];
        next_slope = as[i] + sign * scale * bs[j];
        ++j;
        if (fequals(prev_slope, next_slope))
        {
            breakpoint[num_breaks-1] = next_break;
        }
        else
        {
            breakpoint[num_breaks] = next_break;
            slope[num_breaks] = next_slope;
            prev_slope = next_slope;
            ++num_breaks;
        }
    }
    
    // don't forget the final slope after the last breakpoint
    next_slope = as[i] + sign * scale * bs[j];

    if (fequals(prev_slope, next_slope))
        --num_breaks;
    else
        slope[num_breaks] = next_slope;

    ret->num_breaks = num_breaks;
    ret->intercept = a->intercept + sign * scale * b->intercept;
}


static void
calc_stem_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree)
{
    assert(u != TSK_NULL);
    assert(v != TSK_NULL);
    
    pdata_t *pdata = (pdata_t *)(tree->pdata);

    plf_t *gx = pdata->gx + v;
    plf_t *hx = pdata->hx + v;
    plf_t *gy = pdata->gy + v;
    plf_t *hy = pdata->hy + v;

    double b = 1;
    double minus_b;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    minus_b = -b;

    if (!(tree->node_flags[v] & TSK_NODE_IS_SAMPLE))
    {
        plf_min(gx, b, hx);
        plf_min(gy, b, hy);
    }
    else
    {
        hx->breakpoint[0] = pdata->x[v];
        hx->slope[0] = minus_b;
        hx->slope[1] = b;
        hx->intercept = b*pdata->x[v];
        hx->num_breaks = 1;
        hy->breakpoint[0] = pdata->y[v];
        hy->slope[0] = minus_b;
        hy->slope[1] = b;
        hy->intercept = b*pdata->y[v];
        hy->num_breaks = 1;
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);

    plf_t *gx = pdata->gx + u;
    plf_t *hx = pdata->hx + v;
    plf_t *gy = pdata->gy + u;
    plf_t *hy = pdata->hy + v;
    plf_t *tmp = pdata->tmp;
    
    plf_add(gx, hx, sign, 1, tmp);
    
    
    gx->num_breaks = tmp->num_breaks;
    gx->intercept = tmp->intercept;
    
    memcpy(gx->breakpoint, tmp->breakpoint, tmp->num_breaks*sizeof(double));
    memcpy(gx->slope, tmp->slope, (tmp->num_breaks+1)*sizeof(double));
    
    plf_add(gy, hy, sign, 1, tmp);
    
    gy->num_breaks = tmp->num_breaks;
    gy->intercept = tmp->intercept;

    memcpy(gy->breakpoint, tmp->breakpoint, tmp->num_breaks*sizeof(double));
    memcpy(gy->slope, tmp->slope, (tmp->num_breaks+1)*sizeof(double));
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    if (tree->node_flags[v] & TSK_NODE_IS_SAMPLE)
        return;

    pdata_t *pdata = (pdata_t *)(tree->pdata);

    tsk_id_t u = tree->parent[v];

    plf_t *gv_x = pdata->gx + v;
    plf_t *hv_x = pdata->hx + v;
    plf_t *fv_x = pdata->fx + v;

    plf_t *gv_y = pdata->gy + v;
    plf_t *hv_y = pdata->hy + v;
    plf_t *fv_y = pdata->fy + v;

    if (u == TSK_NULL)
    {
        fv_x->num_breaks = gv_x->num_breaks;
        fv_x->intercept = gv_x->intercept;
        memcpy(fv_x->slope, gv_x->slope,
            (gv_x->num_breaks+1)*sizeof(double));
        memcpy(fv_x->breakpoint, gv_x->breakpoint,
            gv_x->num_breaks*sizeof(double));
        fv_y->num_breaks = gv_y->num_breaks;
        fv_y->intercept = gv_y->intercept;
        memcpy(fv_y->slope, gv_y->slope,
            (gv_y->num_breaks+1)*sizeof(double));
        memcpy(fv_y->breakpoint, gv_y->breakpoint,
            gv_y->num_breaks*sizeof(double));
        return;
    }

    plf_t *fu_x = pdata->fx + u;
    plf_t *fu_y = pdata->fy + u;
    plf_t *tmp = pdata->tmp;

    double b = 1;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    plf_add(fu_x, hv_x, -1, 1, fv_x);
    plf_min(fv_x, b, tmp);
    plf_add(tmp, gv_x, +1, 1, fv_x);

    plf_add(fu_y, hv_y, -1, 1, fv_y);
    plf_min(fv_y, b, tmp);
    plf_add(tmp, gv_y, +1, 1, fv_y);
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
    double len = 0;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        if (!(tree->node_flags[u] & TSK_NODE_IS_SAMPLE))
            len += plf_min_score(pdata->gx + u) + plf_min_score(pdata->gy + u);
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

    plf_t *fx = pdata->fx + node_id;
    plf_t *Fx = mpr_summary->Fx + node_id;
    plf_t *fy = pdata->fy + node_id;
    plf_t *Fy = mpr_summary->Fy + node_id;

    plf_add(fx, Fx, -1, 1, pdata->tmp);
    plf_add(Fx, pdata->tmp, +1, w / total_weight, mpr_summary->tmp);

    Fx->num_breaks = mpr_summary->tmp->num_breaks;
    Fx->intercept = mpr_summary->tmp->intercept;
    memcpy(Fx->slope, mpr_summary->tmp->slope,
        (mpr_summary->tmp->num_breaks+1)*sizeof(double));
    memcpy(Fx->breakpoint, mpr_summary->tmp->breakpoint,
        (mpr_summary->tmp->num_breaks)*sizeof(double));

    plf_add(fy, Fy, -1, 1, pdata->tmp);
    plf_add(Fy, pdata->tmp, +1, w / total_weight, mpr_summary->tmp);

    Fy->num_breaks = mpr_summary->tmp->num_breaks;
    Fy->intercept = mpr_summary->tmp->intercept;
    memcpy(Fy->slope, mpr_summary->tmp->slope,
        (mpr_summary->tmp->num_breaks+1)*sizeof(double));
    memcpy(Fy->breakpoint, mpr_summary->tmp->breakpoint,
        (mpr_summary->tmp->num_breaks)*sizeof(double));
}


SEXP C_treeseq_linear_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP x,
    SEXP y)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int i;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_samples = (int) tsk_treeseq_get_num_samples(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);
    
    plf_t *gx = calloc(num_nodes, sizeof(plf_t));
    plf_t *gy = calloc(num_nodes, sizeof(plf_t));
    plf_t *hx = calloc(num_nodes, sizeof(plf_t));
    plf_t *hy = calloc(num_nodes, sizeof(plf_t));
    plf_t *fx = calloc(num_nodes, sizeof(plf_t));
    plf_t *fy = calloc(num_nodes, sizeof(plf_t));
    plf_t *Fx = calloc(num_nodes, sizeof(plf_t));
    plf_t *Fy = calloc(num_nodes, sizeof(plf_t));

    for (i = 0; i < num_nodes; ++i)
    {
        if (plf_init(gx + i, num_samples)) goto out;
        if (plf_init(gy + i, num_samples)) goto out;
        if (plf_init(hx + i, num_samples)) goto out;
        if (plf_init(hy + i, num_samples)) goto out;
        if (plf_init(fx + i, num_samples)) goto out;
        if (plf_init(fy + i, num_samples)) goto out;
        if (plf_init(Fx + i, num_samples)) goto out;
        if (plf_init(Fy + i, num_samples)) goto out;
    }

    plf_t tmp1;
    plf_t tmp2;

    if (plf_init(&tmp1, num_samples)) goto out;
    if (plf_init(&tmp2, num_samples)) goto out;

    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.Fx = Fx;
    mpr_summary.Fy = Fy;
    mpr_summary.tmp = &tmp2;
    mpr_summary.node_weight = REAL(node_weight);
    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);

    pdata.gx = gx;
    pdata.hx = hx;
    pdata.fx = fx;
    pdata.gy = gy;
    pdata.hy = hy;
    pdata.fy = fy;
    pdata.tmp = &tmp1;
    pdata.x = REAL(x);
    pdata.y = REAL(y);

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

    int k;
    SEXP mpr_x = PROTECT(Rf_allocVector(VECSXP, num_nodes-num_samples));
    SEXP mpr_y = PROTECT(Rf_allocVector(VECSXP, num_nodes-num_samples));

    for (i = num_samples, k = 0; i < num_nodes; ++i, ++k)
    {
        SET_VECTOR_ELT(mpr_x, k, Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_x, k), 0, Rf_ScalarReal(Fx[i].intercept));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_x, k), 1,
            Rf_allocVector(REALSXP, Fx[i].num_breaks+1));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_x, k), 2,
            Rf_allocVector(REALSXP, Fx[i].num_breaks));
        memcpy(REAL(VECTOR_ELT(VECTOR_ELT(mpr_x, k), 1)), Fx[i].slope,
            (Fx[i].num_breaks + 1) * sizeof(double));
        memcpy(REAL(VECTOR_ELT(VECTOR_ELT(mpr_x, k), 2)), Fx[i].breakpoint,
            Fx[i].num_breaks * sizeof(double));

        SET_VECTOR_ELT(mpr_y, k, Rf_allocVector(VECSXP, 3));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_y, k), 0, Rf_ScalarReal(Fy[i].intercept));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_y, k), 1,
            Rf_allocVector(REALSXP, Fy[i].num_breaks+1));
        SET_VECTOR_ELT(VECTOR_ELT(mpr_y, k), 2,
            Rf_allocVector(REALSXP, Fy[i].num_breaks));
        memcpy(REAL(VECTOR_ELT(VECTOR_ELT(mpr_y, k), 1)), Fy[i].slope,
            (Fy[i].num_breaks + 1) * sizeof(double));
        memcpy(REAL(VECTOR_ELT(VECTOR_ELT(mpr_y, k), 2)), Fy[i].breakpoint,
            Fy[i].num_breaks * sizeof(double));
    }

out:
    for (i = 0; i < num_nodes; ++i)
    {
        if (gx) plf_free(gx + i);
        if (gy) plf_free(gy + i);
        if (hx) plf_free(hx + i);
        if (hy) plf_free(hy + i);
        if (hx) plf_free(fx + i);
        if (fy) plf_free(fy + i);
        if (Fx) plf_free(Fx + i);
        if (Fy) plf_free(Fy + i);
    }
    plf_free(&tmp1);
    plf_free(&tmp2);
    free(gx);
    free(gy);
    free(hx);
    free(hy);
    free(fx);
    free(fy);
    free(Fx);
    free(Fy);
    tsx_tree_free(&tree);
    UNPROTECT(4);
    return Rf_list4(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, mpr_x, mpr_y);
}
