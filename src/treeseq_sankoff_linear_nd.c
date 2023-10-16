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
    plf_t **g;
    plf_t **h;
    plf_t **f;
    // workspace
    plf_t *tmp;
    // observed states
    const double *x;
    int num_dims;
} pdata_t;


typedef struct mpr_summary {
    plf_t **F;
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

    plf_t *g;
    plf_t *h;

    int i;
    int num_dims = pdata->num_dims;
    double b = 1;
    double minus_b;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    minus_b = -b;

    if (!(tree->node_flags[v] & TSK_NODE_IS_SAMPLE))
    {
        for (i = 0; i < num_dims; ++i)
        {
            g = pdata->g[i] + v;
            h = pdata->h[i] + v;
            plf_min(g, b, h);
        }
    }
    else
    {
        const double *restrict x = pdata->x + v*num_dims;   
        for (i = 0; i < num_dims; ++i)
        {
            h = pdata->h[i] + v;
            h->breakpoint[0] = x[i];
            h->slope[0] = minus_b;
            h->slope[1] = b;
            h->intercept = b*x[i];
            h->num_breaks = 1;
        }
    }
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);

    int i;
    int num_dims = pdata->num_dims;

    plf_t *g;
    plf_t *h;
    plf_t *tmp = pdata->tmp;
    
    for (i = 0; i < num_dims; ++i)
    {
        g = pdata->g[i] + u;
        h = pdata->h[i] + v;
     
        plf_add(g, h, sign, 1, tmp);

        g->num_breaks = tmp->num_breaks;
        g->intercept = tmp->intercept;
        
        memcpy(g->breakpoint, tmp->breakpoint, tmp->num_breaks*sizeof(double));
        memcpy(g->slope, tmp->slope, (tmp->num_breaks+1)*sizeof(double));
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

    tsk_id_t u = tree->parent[v];

    plf_t *g;
    plf_t *h;
    plf_t *f;

    if (u == TSK_NULL)
    {
        for (i = 0; i < num_dims; ++i)
        {
            g = pdata->g[i] + v;
            f = pdata->f[i] + v;
            f->num_breaks = g->num_breaks;
            f->intercept = g->intercept;
            memcpy(f->slope, g->slope,
                (g->num_breaks+1)*sizeof(double));
            memcpy(f->breakpoint, g->breakpoint,
                g->num_breaks*sizeof(double));
        }
        return;
    }

    plf_t *fu;
    plf_t *tmp = pdata->tmp;

    double b = 1;

    if (tree->time)
        b = 1 / (tree->time[u] - tree->time[v]);

    for (i = 0; i < num_dims; ++i)
    {
        g = pdata->g[i] + v;
        h = pdata->h[i] + v;
        f = pdata->f[i] + v;
        fu = pdata->f[i] + u;
        
        plf_add(fu, h, -1, 1, f);
        plf_min(f, b, tmp);
        plf_add(tmp, g, +1, 1, f);
    }
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
    double len = 0;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        if (!(tree->node_flags[u] & TSK_NODE_IS_SAMPLE))
        {
            for (i = 0; i < num_dims; ++i)
                len += plf_min_score(pdata->g[i] + u);
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

    int i;
    int num_dims = pdata->num_dims;

    double w = t_right - t_left;

    mpr_summary->node_weight[node_id] += w;

    double total_weight = mpr_summary->node_weight[node_id];

    plf_t *f;
    plf_t *F;

    for (i = 0; i < num_dims; ++i)
    {
        f = pdata->f[i] + node_id;
        F = mpr_summary->F[i] + node_id;

        plf_add(f, F, -1, 1, pdata->tmp);
        plf_add(F, pdata->tmp, +1, w / total_weight, mpr_summary->tmp);

        F->num_breaks = mpr_summary->tmp->num_breaks;
        F->intercept = mpr_summary->tmp->intercept;
        memcpy(F->slope, mpr_summary->tmp->slope,
            (mpr_summary->tmp->num_breaks+1)*sizeof(double));
        memcpy(F->breakpoint, mpr_summary->tmp->breakpoint,
            (mpr_summary->tmp->num_breaks)*sizeof(double));

    }
}


SEXP C_treeseq_linear_nd_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP x,
    SEXP nx)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int i;
    int j;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_samples = (int) tsk_treeseq_get_num_samples(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);

    int num_dims = *INTEGER(Rf_getAttrib(x, R_DimSymbol));

    int *num_sites = INTEGER(nx);
    int max_num_sites = 0;
    
    plf_t **g = calloc(num_dims, sizeof(plf_t *));
    plf_t **h = calloc(num_dims, sizeof(plf_t *));
    plf_t **f = calloc(num_dims, sizeof(plf_t *));
    plf_t **F = calloc(num_dims, sizeof(plf_t *));

    
    for (j = 0; j < num_dims; ++j)
    {
        g[j] = calloc(num_nodes, sizeof(plf_t));
        h[j] = calloc(num_nodes, sizeof(plf_t));
        f[j] = calloc(num_nodes, sizeof(plf_t));
        F[j] = calloc(num_nodes, sizeof(plf_t));
        for (i = 0; i < num_nodes; ++i)
        {
            if (plf_init(g[j] + i, num_sites[j])) goto out;
            if (plf_init(h[j] + i, num_sites[j])) goto out;
            if (plf_init(f[j] + i, num_sites[j])) goto out;
            if (plf_init(F[j] + i, num_sites[j])) goto out;
        }
        if (num_sites[j] > max_num_sites)
            max_num_sites = num_sites[j];
    }

    plf_t tmp1;
    plf_t tmp2;

    if (plf_init(&tmp1, max_num_sites)) goto out;
    if (plf_init(&tmp2, max_num_sites)) goto out;

    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    memset(REAL(node_weight), 0, num_nodes * sizeof(double));
    memset(REAL(tree_length), 0, num_trees * sizeof(double));

    mpr_summary.F = F;
    mpr_summary.tmp = &tmp2;
    mpr_summary.node_weight = REAL(node_weight);
    mpr_summary.tree_weight = 0;
    mpr_summary.mean_tree_length = 0;
    mpr_summary.tree_length = REAL(tree_length);

    pdata.g = g;
    pdata.h = h;
    pdata.f = f;
    pdata.tmp = &tmp1;
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

    int k;
    SEXP mpr = PROTECT(Rf_allocVector(VECSXP, num_dims));

    for (j = 0; j < num_dims; ++j)
    {
        SET_VECTOR_ELT(
            mpr
            , j
            , Rf_allocVector(VECSXP, num_nodes-num_samples)
        );
        for (i = num_samples, k = 0; i < num_nodes; ++i, ++k)
        {
            SET_VECTOR_ELT(
                VECTOR_ELT(mpr, j)
                , k
                , Rf_allocVector(VECSXP, 3)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(mpr, j), k)
                , 0
                , Rf_ScalarReal((F[j] + i)->intercept)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(mpr, j), k)
                , 1
                , Rf_allocVector(REALSXP, (F[j] + i)->num_breaks + 1)
            );
            SET_VECTOR_ELT(
                VECTOR_ELT(VECTOR_ELT(mpr, j), k)
                , 2
                , Rf_allocVector(REALSXP, (F[j] + i)->num_breaks)
            );
            memcpy(
                REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(mpr, j), k), 1))
                , (F[j] + i)->slope
                , ((F[j] + i)->num_breaks + 1) * sizeof(double)
            );
            memcpy(
                REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(mpr, j), k), 2))
                , (F[j] + i)->breakpoint
                , (F[j] + i)->num_breaks * sizeof(double)
            );
        }
    }

out:
    
    for (j = 0; j < num_dims; ++j)
    {
        for (i = 0; i < num_nodes; ++i)
        {
            if (g[j]) plf_free(g[j] + i);
            if (h[j]) plf_free(h[j] + i);
            if (f[j]) plf_free(f[j] + i);
            if (F[j]) plf_free(F[j] + i);
        }
        free(g[j]);
        free(h[j]);
        free(f[j]);
        free(F[j]);
    }
    plf_free(&tmp1);
    plf_free(&tmp2);
    free(g);
    free(h);
    free(f);
    free(F);
    tsx_tree_free(&tree);
    UNPROTECT(3);
    return Rf_list3(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, mpr);
}


SEXP C_treeseq_linear_nd_mpr_minimize(SEXP F)
{
    int i;
    int j;
    int k;
    int l;
    int num_dims = Rf_length(F);
    int n = Rf_length(VECTOR_ELT(F, 0));

    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, num_dims));
    double *restrict x = REAL(ans);
    double *restrict a;
    double *restrict b;
    GetRNGstate();
    for (j = 0; j < num_dims; ++j)
    {
        for (i = 0; i < n; ++i)
        {
            a = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 1));
            b = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
            k = 0;
            while (a[k+1] < 0)
                ++k;
            x[i+j*n] = b[k];
            l = 1;
            while (a[k+1] == 0)
            {
                ++l;
                ++k;
                if (unif_rand() < (1 / (double)l))
                    x[i+j*n] = b[k];
            }
        }
    }
    PutRNGstate();
    UNPROTECT(1);
    return ans;
}


static double
compute_score(double x, double intercept, double *slope, double *breakpoint,
    int num_breaks)
{
    int j;
    double score = 0;
    if (x < breakpoint[0])
    {
        score = intercept + x*slope[0];
    }
    else
    {
        j = 0;
        score = intercept + slope[0]*breakpoint[0];
        while (j < (num_breaks-1) && x > breakpoint[j+1])
        {
            ++j;
            score += slope[j] * (breakpoint[j] - breakpoint[j-1]);
        }
        score += slope[j+1] * (x - breakpoint[j]);
    }
    return score;
}


SEXP C_treeseq_linear_nd_mpr_minimize_discrete(SEXP F, SEXP sites)
{
    int i;
    int j;
    int k;
    int num_dims = Rf_length(F);
    int n = Rf_length(VECTOR_ELT(F, 0));
    int num_breaks;
    int num_sites = *INTEGER(Rf_getAttrib(sites, R_DimSymbol));
    
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, num_dims));
    
    double *restrict x = REAL(ans);
    double *restrict coords = REAL(sites);

    double intercept;
    double *restrict a;
    double *restrict b;

    double score;
    double min_score;

    for (i = 0; i < n; ++i)
    {
        min_score = R_PosInf;
        for (k = 0; k < num_sites; ++k)
        {
            score = 0;
            for (j = 0; j < num_dims; ++j)
            {
                intercept = *REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 0));
                a = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 1));
                b = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
                num_breaks = Rf_length(
                    VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(F, j), i), 2));
                score += compute_score(
                    coords[k+j*num_sites], intercept, a, b, num_breaks);
            }
            if (score < min_score)
            {
                min_score = score;
                for (j = 0; j < num_dims; ++j)
                    x[i+j*n] = coords[k+j*num_sites];
            }
        }
    }
    UNPROTECT(1);
    return ans;
}
