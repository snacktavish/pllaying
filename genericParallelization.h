#ifndef _GENERIC_PARALL_H 
#define _GENERIC_PARALL_H 


/**********/
/* CONFIG */
/**********/
#define DEBUG_PARALLEL
#define GENERIC_PARALLELIZATION

#define NOT ! 
#define IS_PARALLEL (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI)) 
void *likelihoodThread(void *tData); 



/******************/
/* MPI SPECIFIC   */
/******************/
#ifdef _FINE_GRAIN_MPI
#include <mpi.h>
#include "mpiHelpers.h"
#define MASTER_P (processID == 0)
#define POP_OR_PUT(buf, elem)  (MASTER_P ? (addIntToBuf((buf++), (elem)) , printf("SEND %d\n", *elem ) ) : ( popIntFromBuf((buf++), (elem)), printf("RECV %d\n", *elem))) 

#define ASSIGN_INT(x,y) (MPI_Bcast(&y,1,MPI_INT,0,MPI_COMM_WORLD),printf("SEND/RECV %d\n", y)); 
#define ASSIGN_BUF(x,y) ((POP_OR_PUT(bufPtr, &y)),assertCtr++)
#define ASSIGN_DBL(x,y) (MPI_Bcast(&y,1,MPI_DOUBLE, 0, MPI_COMM_WORLD), printf("SEND/RECV %d\n", y)); 

extern int processes; 
extern int processID; 

int* addIntToBuf(int* buf, int *toAdd);
int* popIntFromBuf(int *buf, int *result); 
#endif 

/*********************/
/* PTHREAD SPECIFIC  */
/*********************/
#ifdef _USE_PTHREADS
#include <pthread.h>
#define MASTER_P (tid == 0)
#define ASSIGN_INT(x,y) (x = y)
#define ASSIGN_BUF(x,y) (x = y)
#define ASSIGN_DBL(x,y) (x = y)
#endif

#endif	/* end include guard  */
