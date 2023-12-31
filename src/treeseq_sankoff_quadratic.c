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
**   F(x,y) = p1*x*x + p2*y*y + p3*x + p4*y + p5
**
** Which gives the minimum sum of squared distances between all
** ancestor-descendant pairs needed to explain the spatial
** distribution of samples assuming the non-sample node is in state
** (x, y).
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
    const double *y;
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
    int ofs = v*5;
    double b = 1;
    double *restrict h_v = pdata->h + ofs;

    if (tree->time)
    {
        b = tree->time[u] - tree->time[v];
    }
    
    if (!(tree->node_flags[v] & TSK_NODE_IS_SAMPLE))
    {
        // internal node
        double *restrict g_v = pdata->g + ofs;
        double p1 = g_v[0];
        double p2 = g_v[1];
        double p3 = g_v[2];
        double p4 = g_v[3];
        double p5 = g_v[4];

        h_v[0] = p1 / (b*p1 + 1);
        h_v[1] = p2 / (b*p2 + 1);
        h_v[2] = p3 / (b*p1 + 1);
        h_v[3] = p4 / (b*p2 + 1);
        h_v[4] = - ( p3*p3 ) / (4 * (p1 + 1/b) )
                 - ( p4*p4 ) / (4 * (p2 + 1/b) ) 
                 + p5;
    }
    else
    {
        // terminal node
        double x = pdata->x[v];
        double y = pdata->y[v];

        h_v[0] = 1 / b;
        h_v[1] = 1 / b;
        h_v[2] = (-2*x) / b;
        h_v[3] = (-2*y) / b;
        h_v[4] = (x*x + y*y) / b;
    }
}


static void
calc_final_cost(tsk_id_t v, tsx_tree_t *tree)
{
    if (tree->node_flags[v] & TSK_NODE_IS_SAMPLE)
        return;

    pdata_t *pdata = (pdata_t *)(tree->pdata);
    int ofs = v*5;
    tsk_id_t u = tree->parent[v];
    if (u == TSK_NULL)
    {
        memcpy(pdata->f + ofs, pdata->g + ofs, 5*sizeof(double));
        return;
    }
    double *restrict g_v = pdata->g + ofs;
    double *restrict h_v = pdata->h + ofs;
    double *restrict f_v = pdata->f + ofs;
    double *restrict f_u = pdata->f + u*5;

    double p1 = f_u[0] - h_v[0];
    double p2 = f_u[1] - h_v[1];
    double p3 = f_u[2] - h_v[2];
    double p4 = f_u[3] - h_v[3];
    double p5 = f_u[4] - h_v[4];

    double b = 1;

    if (tree->time)
    {
        b = tree->time[u] - tree->time[v];
    }

    f_v[0] = p1 / (b*p1 + 1);
    f_v[1] = p2 / (b*p2 + 1);
    f_v[2] = p3 / (b*p1 + 1);
    f_v[3] = p4 / (b*p2 + 1);
    f_v[4] = - ( p3*p3 ) / (4 * (p1 + 1/b) )
             - ( p4*p4 ) / (4 * (p2 + 1/b) )
             + p5;

    f_v[0] += g_v[0];
    f_v[1] += g_v[1];
    f_v[2] += g_v[2];
    f_v[3] += g_v[3];
    f_v[4] += g_v[4];
}


static void
increment_node_cost(tsk_id_t u, tsk_id_t v, tsx_tree_t *tree, int sign)
{
    pdata_t *pdata = (pdata_t *)(tree->pdata);
    double *restrict g_u = pdata->g + u*5;
    double *restrict h_v = pdata->h + v*5;
    g_u[0] += sign*h_v[0];
    g_u[1] += sign*h_v[1];
    g_u[2] += sign*h_v[2];
    g_u[3] += sign*h_v[3];
    g_u[4] += sign*h_v[4];
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
    double *g;
    double p1, p2, p3, p4, p5, min_s, len = 0;
    for (u = left_child[virtual_root]; u != TSK_NULL; u = right_sib[u])
    {
        if (!(tree->node_flags[u] & TSK_NODE_IS_SAMPLE))
        {
            g = pdata->g + u*5;
            p1 = g[0];
            p2 = g[1];
            p3 = g[2];
            p4 = g[3];
            p5 = g[4];
            min_s = p5 - (p3*p3)/(4*p1) - (p4*p4)/(4*p2);
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

    const double *restrict f = pdata->f + 5 * node_id;
    double *restrict F = mpr_summary->F + 5 * node_id;

    F[0] += w * (f[0] - F[0]) / total_weight;
    F[1] += w * (f[1] - F[1]) / total_weight;
    F[2] += w * (f[2] - F[2]) / total_weight;
    F[3] += w * (f[3] - F[3]) / total_weight;
    F[4] += w * (f[4] - F[4]) / total_weight;
}


SEXP C_treeseq_quadratic_mpr(
    SEXP treeseq,
    SEXP use_brlen, 
    SEXP x,
    SEXP y)
{
    tsx_tree_t tree;
    pdata_t pdata;
    mpr_summary_t mpr_summary;
    tsk_treeseq_t *ts = (tsk_treeseq_t *) R_ExternalPtrAddr(treeseq);
    int nprotect = 0;
    int num_nodes = (int) tsk_treeseq_get_num_nodes(ts);
    int num_trees = (int) tsk_treeseq_get_num_trees(ts);
    
    SEXP g = PROTECT(Rf_allocVector(REALSXP, 5*num_nodes));
    SEXP h = PROTECT(Rf_allocVector(REALSXP, 5*num_nodes));
    SEXP f = PROTECT(Rf_allocVector(REALSXP, 5*num_nodes));

    SEXP F = PROTECT(Rf_allocMatrix(REALSXP, 5, num_nodes));
    SEXP node_weight = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP tree_length = PROTECT(Rf_allocVector(REALSXP, num_trees));

    nprotect = 6;

    memset(REAL(g), 0, 5 * num_nodes * sizeof(double));
    memset(REAL(h), 0, 5 * num_nodes * sizeof(double));
    memset(REAL(f), 0, 5 * num_nodes * sizeof(double));
    memset(REAL(F), 0, 5 * num_nodes * sizeof(double));
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

out:
    tsx_tree_free(&tree);
    UNPROTECT(nprotect);
    return Rf_list3(
        Rf_ScalarReal(mpr_summary.mean_tree_length), tree_length, F);
}


SEXP C_treeseq_quadratic_mpr_minimize(SEXP F)
{
    int i;
    int nr = *INTEGER(Rf_getAttrib(F, R_DimSymbol));
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, nr, 2));
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0,
        VECTOR_ELT(Rf_getAttrib(F, R_DimNamesSymbol), 0));
    SET_VECTOR_ELT(dimnames, 1, R_NilValue);
    double *x = REAL(ans);
    double *y = REAL(ans) + nr;
    const double *restrict p = REAL(F);
    for (i = 0; i < nr; ++i)
    {
        x[i] = -p[i+2*nr] / (2*p[i+0*nr]);
        y[i] = -p[i+3*nr] / (2*p[i+1*nr]);
    }
    Rf_setAttrib(ans, R_DimNamesSymbol, dimnames);
    UNPROTECT(2);
    return ans;
}


SEXP C_treeseq_quadratic_mpr_sample(SEXP F, SEXP rate)
{
    int i;
    int sign;
    int nr = *INTEGER(Rf_getAttrib(F, R_DimSymbol));

    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, nr, 2));
    SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0,
        VECTOR_ELT(Rf_getAttrib(F, R_DimNamesSymbol), 0));
    SET_VECTOR_ELT(dimnames, 1, R_NilValue);
    
    double *x = REAL(ans);
    double *y = REAL(ans) + nr;
    const double *restrict p = REAL(F);

    double lambda = *REAL(rate);

    double r;
    double phi;
    double cosphi;
    double sinphi;
    double xp;
    double yp;
    double min_x;
    double min_y;
    double min_score;
    double score;
    double step;
    double step_x;
    double step_y;
    GetRNGstate();
    for (i = 0; i < nr; ++i)
    {
        min_x = -p[i+2*nr] / (2*p[i+0*nr]);
        min_y = -p[i+3*nr] / (2*p[i+1*nr]);
        min_score = p[i+4*nr] - 
                    (p[i+2*nr]*p[i+2*nr]) / (4*p[i+0*nr]) - 
                    (p[i+3*nr]*p[i+3*nr]) / (4*p[i+1*nr]);
        xp = min_x;
        yp = min_y;
        step = 1;
        sign = +1;

        r = -log(1 - unif_rand()) / lambda; // rexp(lambda)
        
        phi = unif_rand() * M_2PI;
        cosphi = cos(phi);
        sinphi = sin(phi);
        step_x = cosphi;
        step_y = sinphi;

        // we want to sample from exp(-(score - min_score)). i.e., from
        // a rate 1 exponential distribution. starting from the minimum,
        // search along a line until (score - min_score) ≈ 5, which is the
        // distance that will include about 99 percent of the exponential 
        // distribution.
        do {
            xp += sign*step_x;
            yp += sign*step_y;
            score = p[i+0*nr]*xp*xp + 
                    p[i+1*nr]*yp*yp + 
                    p[i+2*nr]*xp + 
                    p[i+3*nr]*yp + 
                    p[i+4*nr];
            if ((score - min_score) > r)
            {
                sign = -1;
                step /= 2;
                step_x = step * cosphi;
                step_y = step * sinphi;
            }
            else if ((score - min_score) < r)
            {
                sign = +1;
                step *= 2;
                step_x = step * cosphi;
                step_y = step * sinphi;
            }
        } while (fabs((score - min_score) - r) > 1e-8);

        x[i] = xp;
        y[i] = yp;

    }
    PutRNGstate();
    Rf_setAttrib(ans, R_DimNamesSymbol, dimnames);
    UNPROTECT(2);
    return ans;
}


SEXP C_treeseq_quadratic_mpr_minimize_discrete(
    SEXP F, SEXP eastings, SEXP northings, SEXP raster)
{
    int i;
    int j;
    int k;
    int l;
    int n = *INTEGER(Rf_getAttrib(F, R_DimSymbol));
    
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, 2));
    double *x = REAL(ans);
    double *y = REAL(ans) + n;

    int start;
    int stop;
    int *ri = INTEGER(R_do_slot(raster, Rf_mkString("i")));
    int *rj = INTEGER(R_do_slot(raster, Rf_mkString("p")));

    const double *restrict p = REAL(F);

    int ncol = Rf_length(eastings);

    double score;
    double min_score;
    double min_x;
    double min_y;

    double *easting = REAL(eastings);
    double *northing = REAL(northings);

    for (l = 0; l < n; ++l)
    {
        min_score = R_PosInf;

        // loop over columns
        for (j = 0; j < ncol; ++j)
        {
            start = rj[j];
            stop = start + rj[j+1] - rj[j];
            for (k = start; k < stop; ++k)
            {
                i = ri[k];
                
                // compute the score for (i,j) cell
                score = p[l+0*n]*easting[j]*easting[j] + 
                    p[l+1*n]*northing[i]*northing[i] + 
                    p[l+2*n]*easting[j] + 
                    p[l+3*n]*northing[i] + 
                    p[l+4*n];

                if (score < min_score)
                {
                    min_x = easting[j];
                    min_y = northing[i];
                    min_score = score;
                }
            }
        }
        x[l] = min_x;
        y[l] = min_y;
    }

    UNPROTECT(1);
    return ans;
}


SEXP C_treeseq_quadratic_mpr_sample_discrete(
    SEXP F, SEXP eastings, SEXP northings, SEXP raster)
{
    int i;
    int j;
    int k;
    int l;
    int n = *INTEGER(Rf_getAttrib(F, R_DimSymbol));
    
    SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n, 2));
    double *x = REAL(ans);
    double *y = REAL(ans) + n;

    int start;
    int stop;
    int *ri = INTEGER(R_do_slot(raster, Rf_mkString("i")));
    int *rj = INTEGER(R_do_slot(raster, Rf_mkString("p")));

    const double *restrict p = REAL(F);

    int ncol = Rf_length(eastings);

    double score;
    double max_score;
    double gumbel;
    double min_x;
    double min_y;

    double *easting = REAL(eastings);
    double *northing = REAL(northings);

    GetRNGstate();
    for (l = 0; l < n; ++l)
    {
        max_score = R_NegInf;

        // loop over columns
        for (j = 0; j < ncol; ++j)
        {
            start = rj[j];
            stop = start + rj[j+1] - rj[j];
            for (k = start; k < stop; ++k)
            {
                i = ri[k];
                
                gumbel = -log(-log(unif_rand()));
                
                // compute the score for (i,j) cell
                score = p[l+0*n]*easting[j]*easting[j] + 
                    p[l+1*n]*northing[i]*northing[i] + 
                    p[l+2*n]*easting[j] + 
                    p[l+3*n]*northing[i] + 
                    p[l+4*n];

                if ((gumbel - score) > max_score)
                {
                    min_x = easting[j];
                    min_y = northing[i];
                    max_score = gumbel - score;
                }
            }
        }
        x[l] = min_x;
        y[l] = min_y;
    }
    PutRNGstate();

    UNPROTECT(1);
    return ans;
}