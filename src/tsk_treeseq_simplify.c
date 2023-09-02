#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <tskit.h>
#include "error.h"

// defined in tsx_treeseq_load.c
void
tsx_treeseq_exptr_free(SEXP treeseq);

void
tsx_tree_exptr_free(SEXP treeseq);

SEXP
C_treeseq_simplify(SEXP treeseq, SEXP sample_ids)
{
    int ret;
    tsk_size_t num_samples = (tsk_size_t)Rf_length(sample_ids);
    tsk_treeseq_t *ts = (tsk_treeseq_t *)R_ExternalPtrAddr(treeseq);
    tsk_size_t N = tsk_treeseq_get_num_nodes(ts);
    const tsk_id_t *samples = (tsk_id_t *)INTEGER(sample_ids);
    tsk_flags_t options = TSK_SIMPLIFY_FILTER_SITES |
                          TSK_SIMPLIFY_FILTER_POPULATIONS |
                          TSK_SIMPLIFY_FILTER_INDIVIDUALS |
                          TSK_SIMPLIFY_REDUCE_TO_SITE_TOPOLOGY;
    tsk_treeseq_t *tss = malloc(sizeof(*tss));
    tsk_tree_t *tree = malloc(sizeof(*tree));
    tsk_id_t *node_map = malloc(N*sizeof(tsk_id_t));
    if (!tss || !tree || !node_map)
        TSX_ERROR("memory allocation failed");
    ret = tsk_treeseq_simplify(ts,samples,num_samples,options,tss,node_map);
    check_tsk_error(ret);
    ret = tsk_tree_init(tree, tss, 0);
    check_tsk_error(ret);
    SEXP exptr1 = PROTECT(R_MakeExternalPtr(tss, NULL, NULL));
    R_RegisterCFinalizer(exptr1, tsx_treeseq_exptr_free);
    SEXP exptr2 = PROTECT(R_MakeExternalPtr(tree, NULL, NULL));
    R_RegisterCFinalizer(exptr2, tsx_tree_exptr_free);
    SEXP nodemap = PROTECT(Rf_allocVector(INTSXP, N));
    memcpy(INTEGER(nodemap), (int *)node_map, N * sizeof(int));
    UNPROTECT(3);
    free(node_map);
    return Rf_list3(exptr1, exptr2, nodemap);
}
