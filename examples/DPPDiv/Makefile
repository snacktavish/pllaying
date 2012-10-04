CPP = g++
CXXFLAGS = -DPHYLO_LIBRARY -c -O2
CC  = gcc
DPPDIV_SRC_DIR = dppdiv_src/
PHYLOLIB_SRC_DIR = ../../
OBJ_DIR = obj/
CXXFILES = Alignment.cpp Calibration.cpp dppdiv.cpp MbEigensystem.cpp MbMath.cpp MbRandom.cpp MbTransitionMatrix.cpp Mcmc.cpp Model.cpp Parameter_basefreq.cpp Parameter_cphyperp.cpp Parameter.cpp Parameter_exchangeability.cpp Parameter_expcalib.cpp Parameter_rate.cpp Parameter_shape.cpp Parameter_speciaton.cpp Parameter_tree.cpp Parameter_treescale.cpp
LIB_CFILES   = bipartitionList.c evaluateGenericSpecial.c evaluatePartialGenericSpecial.c makenewzGenericSpecial.c models.c optimizeModel.c restartHashTable.c searchAlgo.c topologies.c trash.c treeIO.c newviewGenericSpecial.c fastDNAparsimony.c randomTree.c recom.c utils.c
PHYLOLIB_OBJ = $(patsubst %.c, $(PHYLOLIB_SRC_DIR)%.o, $(LIB_CFILES))
DPPDIV_OBJ   = $(patsubst %.cpp, $(DPPDIV_SRC_DIR)%.o, $(CXXFILES))
CFLAGS       = -I$(PHYLOLIB_SRC_DIR) -c -g -march=native -O2 -fomit-frame-pointer -funroll-loops -D_GNU_SOURCE -D__SIM_SSE3 -D_OPTIMIZED_FUNCTIONS -msse3 -g

all: phylolib dppdiv

phylolib: $(PHYLOLIB_OBJ)
	mv $(PHYLOLIB_SRC_DIR)*.o $(OBJ_DIR)

dppdiv: $(DPPDIV_OBJ)
	mv $(DPPDIV_SRC_DIR)*.o $(OBJ_DIR)

%.o: %.cpp
	$(CPP) $(CXXFLAGS)  $+ -I$(DPPDIV_SRC_DIR) -o $@

%.o: $(PHYLOLIB_SRC_DIR)%.c
	$(CC) $(CFLAGS) -c $+ -o $@

