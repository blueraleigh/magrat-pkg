#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <tskit.h>

#include "error.h"


/*
** Given a flat index (k) into the lower triangle (excluding the diagonal) 
** of a n-by-n matrix stored in column major order, return the corresponding 
** row and column indices (i,j).
**
** Modified from Algorithm 3 of the following reference
**   https://hal.archives-ouvertes.fr/hal-02047514/document
**
** The linked reference uses a 1-based index. 
** Here we use a 0-based index.
**
** See also: https://en.wikipedia.org/wiki/Triangular_number
*/
static void
k2ij(int n, int k, tsk_id_t *i, tsk_id_t *j)
{
    ++k;
    --n;
    int kp = ( n * (n + 1) ) / 2 - k;
    int p = ( sqrt(1 + 8 * kp) - 1 ) / 2;
    *i = k - n*(n-1)/2 + p*(p+1)/2;
    *j = n - p - 1;
}


static void
tsx_treeseq_tmrca(tsk_tree_t *tree, double *restrict tmrca)
{
    int ret;
    int k;
    tsk_id_t i;
    tsk_id_t j;
    tsk_id_t mrca;
    tsk_id_t num_samples = (tsk_id_t)(tree->tree_sequence->num_samples);
    double tree_span;
    const double *restrict time = tree->tree_sequence->tables->nodes.time;

    int n = ( num_samples * (num_samples  - 1) ) / 2;

    double *restrict weight = calloc(n, sizeof(*weight));

    if (!weight)
    {
        TSX_WARN("memory allocation failed");
        goto out;
    }

    for (
        ret = tsk_tree_first(tree);
        ret == TSK_TREE_OK;
        ret = tsk_tree_next(tree))
    {
        tree_span = tree->interval.right - tree->interval.left;
        for (k = 0; k < n; ++k)
        {
            k2ij(num_samples, k, &i, &j);
            tsk_tree_get_mrca(tree, i, j, &mrca);
            if (mrca != TSK_NULL)
            {
                weight[k] += tree_span;
                tmrca[k] += tree_span * (time[mrca] - tmrca[k]) / weight[k];
            }
        }
    }
    for (k = 0; k < n; ++k)
    {
        if (weight[k] == 0)
        {
            tmrca[k] = NA_REAL;
        }
    }
out:
    free(weight);
}


SEXP C_treeseq_tmrca(SEXP tr)
{
    tsk_tree_t *tree = (tsk_tree_t *)R_ExternalPtrAddr(tr);
    int num_samples = (int)(tree->tree_sequence->num_samples);
    int n = (num_samples * (num_samples - 1)) / 2;
    SEXP ret = PROTECT(Rf_allocVector(REALSXP, n));
    double *tmrca = REAL(ret);
    memset(tmrca, 0, n * sizeof(*tmrca));
    tsx_treeseq_tmrca(tree, tmrca);
    UNPROTECT(1);
    return ret;
}

