CC = @MPICC@

lib_LTLIBRARIES =
lib_LIBRARIES =

if BUILD_MPI 

if BUILD_SSE3
lib_LTLIBRARIES += libpll-sse3-mpi.la
lib_LIBRARIES += libpll-sse3-mpi.a
libpll_sse3_mpi_la_SOURCES = $(top_srcdir)/src/genericParallelization.c $(top_srcdir)/src/hash.c $(top_srcdir)/src/stack.c $(top_srcdir)/src/ssort.c $(top_srcdir)/src/queue.c $(top_srcdir)/src/utils.c $(top_srcdir)/src/randomTree.c $(top_srcdir)/src/optimizeModel.c $(top_srcdir)/src/trash.c $(top_srcdir)/src/searchAlgo.c $(top_srcdir)/src/topologies.c $(top_srcdir)/src/fastDNAparsimony.c $(top_srcdir)/src/treeIO.c $(top_srcdir)/src/models.c $(top_srcdir)/src/evaluatePartialGenericSpecial.c $(top_srcdir)/src/evaluateGenericSpecial.c $(top_srcdir)/src/newviewGenericSpecial.c $(top_srcdir)/src/makenewzGenericSpecial.c $(top_srcdir)/src/bipartitionList.c $(top_srcdir)/src/restartHashTable.c $(top_srcdir)/src/recom.c $(top_srcdir)/src/lexer.c $(top_srcdir)/src/alignment.c $(top_srcdir)/src/newick.c $(top_srcdir)/src/parsePartition.c
libpll_sse3_mpi_la_CFLAGS = -D_FINE_GRAIN_MPI -D__SSE3 -msse3 -c -O2 -fomit-frame-pointer -funroll-loops -Wall -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef 
libpll_sse3_mpi_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
libpll_sse3_mpi_a_SOURCES = $(libpll_sse3_mpi_la_SOURCES)
libpll_sse3_mpi_a_CFLAGS = $(libpll_sse3_mpi_la_CFLAGS)
endif

if BUILD_AVX
lib_LTLIBRARIES += libpll-avx-mpi.la
lib_LIBRARIES += libpll-avx-mpi.a
libpll_avx_mpi_la_SOURCES = $(top_srcdir)/src/avxLikelihood.c $(top_srcdir)/src/genericParallelization.c $(top_srcdir)/src/hash.c $(top_srcdir)/src/stack.c $(top_srcdir)/src/ssort.c $(top_srcdir)/src/queue.c $(top_srcdir)/src/utils.c $(top_srcdir)/src/randomTree.c $(top_srcdir)/src/optimizeModel.c $(top_srcdir)/src/trash.c $(top_srcdir)/src/searchAlgo.c $(top_srcdir)/src/topologies.c $(top_srcdir)/src/fastDNAparsimony.c $(top_srcdir)/src/treeIO.c $(top_srcdir)/src/models.c $(top_srcdir)/src/evaluatePartialGenericSpecial.c $(top_srcdir)/src/evaluateGenericSpecial.c $(top_srcdir)/src/newviewGenericSpecial.c $(top_srcdir)/src/makenewzGenericSpecial.c $(top_srcdir)/src/bipartitionList.c $(top_srcdir)/src/restartHashTable.c $(top_srcdir)/src/recom.c $(top_srcdir)/src/lexer.c $(top_srcdir)/src/alignment.c $(top_srcdir)/src/newick.c $(top_srcdir)/src/parsePartition.c
libpll_avx_mpi_la_CFLAGS = -D_FINE_GRAIN_MPI -D__AVX -mavx -D__SSE3 -msse3 -c -O2 -fomit-frame-pointer -funroll-loops -Wall -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes -Wpointer-sign -Wextra -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable -Wformat  -Wformat-nonliteral -Wparentheses -Wsequence-point -Wuninitialized -Wundef 
libpll_avx_mpi_la_LDFLAGS = -version-info @LIBPLL_MAJOR@:@LIBPLL_MINOR@:@LIBPLL_REV@
libpll_avx_mpi_a_SOURCES = $(libpll_avx_mpi_la_SOURCES)
libpll_avx_mpi_a_CFLAGS = $(libpll_avx_mpi_la_CFLAGS)
endif


endif
