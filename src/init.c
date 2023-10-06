#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define FDEF(name, n)  {#name, (DL_FUNC) &F77_SYMBOL(name), n, NULL}

void
F77_NAME(spheremean)(int *,double *,double *,double *,double *,double *,
    double *, double *);

static R_FortranMethodDef fortranMethods[] = {
    {"F_spheremean", (DL_FUNC) &F77_NAME(spheremean), 8},

    {NULL, NULL, 0}
};


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
SEXP C_treeseq_sankoff_island_mpr_sample(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_treeseq_sankoff_planar_mpr(SEXP,SEXP,SEXP,SEXP);

SEXP C_treeseq_tmrca(SEXP);

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
    CALLDEF(C_treeseq_sankoff_island_mpr_sample, 5),
    CALLDEF(C_treeseq_sankoff_planar_mpr, 4),

    CALLDEF(C_treeseq_tmrca, 1),

    {NULL, NULL, 0}
};


void attribute_visible R_init_magrat(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, fortranMethods, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
