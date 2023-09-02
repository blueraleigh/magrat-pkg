CC = cc

CFLAGS = -Wall -std=c99 -pedantic -fPIC

CFLAGS += -I.

CFLAGS += -I./subprojects/kastore

TSK_OBJECTS = kastore.o tskit_convert.o tskit_core.o tskit_genotypes.o \
	tskit_haplotype_matching.o tskit_stats.o tskit_tables.o tskit_trees.o 

libtskit.a: ${TSK_OBJECTS}
	${AR} rcs $@ ${TSK_OBJECTS}

kastore.o:
	${CC} -c ${CFLAGS} ./subprojects/kastore/kastore.c -o kastore.o

tskit_%.o:
	${CC} -c ${CFLAGS} ./tskit/$*.c -o $@
