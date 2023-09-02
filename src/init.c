#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP C_treeseq_load(SEXP);
SEXP C_treeseq_to_phylo(SEXP);
SEXP C_treeseq_simplify(SEXP,SEXP);
SEXP C_treeseq_sample(SEXP);
SEXP C_treeseq_nodes(SEXP);
SEXP C_treeseq_edges(SEXP);
SEXP C_treeseq_individuals(SEXP);
SEXP C_treeseq_intervals(SEXP);

SEXP C_treeseq_sankoff_lattice_mpr(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_sankoff_island_mpr(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_sankoff_planar_mpr(SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_treeseq_load, 1),
    CALLDEF(C_treeseq_to_phylo, 1),
    CALLDEF(C_treeseq_simplify, 2),
    CALLDEF(C_treeseq_sample, 1),
    CALLDEF(C_treeseq_nodes, 1),
    CALLDEF(C_treeseq_edges, 1),
    CALLDEF(C_treeseq_individuals, 1),
    CALLDEF(C_treeseq_intervals, 1),

    CALLDEF(C_treeseq_sankoff_lattice_mpr, 8),
    CALLDEF(C_treeseq_sankoff_island_mpr, 5),
    CALLDEF(C_treeseq_sankoff_planar_mpr, 4),

    {NULL, NULL, 0}
};


void attribute_visible R_init_magrat(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
