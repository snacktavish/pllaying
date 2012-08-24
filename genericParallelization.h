#ifndef _GENERIC_PARALL_H 
#define _GENERIC_PARALL_H 


extern double *globalResult; 


/**********/
/* CONFIG */
/**********/

#define _PORTABLE_PTHREADS
/* #define DEBUG_PARALLEL */
/* #define DEBUG_MPI_EACH_SEND */



#define NOT ! 
#define IS_PARALLEL (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI)) 
void *likelihoodThread(void *tData); 



/******************/
/* MPI SPECIFIC   */
/******************/
#ifdef _FINE_GRAIN_MPI
#include <mpi.h>
#include "mpiHelpers.h"

#ifdef DEBUG_MPI_EACH_SEND
#define DEBUG_PRINT(text, elem) printf(text, elem)
#else 
#define DEBUG_PRINT(text, elem) NULL
#endif

#define VOLATILE_PAR 
#define MASTER_P (processID == 0)
#define POP_OR_PUT(buf, elem)  (MASTER_P ? (addIntToBuf((buf++), (elem)) , DEBUG_PRINT("\tSEND %d\n", *elem ) ) : ( popIntFromBuf((buf++), (elem)), DEBUG_PRINT("\tRECV %d\n", *elem))) 
#define POP_OR_PUT_DBL(buf, elem) (MASTER_P ? (addDblToBuf((buf++), (elem)) , DEBUG_PRINT("\tSEND %f\n", *elem ) ) : ( popDblFromBuf((buf++), (elem)), DEBUG_PRINT("\tRECV %f\n", *elem))) 

#define ASSIGN_INT(x,y) (MPI_Bcast(&y,1,MPI_INT,0,MPI_COMM_WORLD),DEBUG_PRINT("\tSEND/RECV %d\n", y)) 
#define ASSIGN_BUF(x,y) ((POP_OR_PUT(bufPtr, &y)),assertCtr++)
#define ASSIGN_BUF_DBL(x,y) (POP_OR_PUT_DBL(bufPtrDbl,&y))
#define ASSIGN_DBL(x,y) (MPI_Bcast(&y,1,MPI_DOUBLE, 0, MPI_COMM_WORLD), DEBUG_PRINT("\tSEND/RECV %f\n", y)) 
#define ASSIGN_DBLS(tar,src,length) MPI_Bcast(tar, length, MPI_DOUBLE, 0, MPI_COMM_WORLD)
#define DOUBLE MPI_DOUBLE
#define ASSIGN_GATHER(tar,src,length,type,tid) MPI_Gather(src,length,type,tar,length,type,0, MPI_COMM_WORLD)
#define SEND_BUF(buf, bufSize,type) if(MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define RECV_BUF(buf, bufSize,type) if(NOT MASTER_P) MPI_Bcast(buf, bufSize, type, 0, MPI_COMM_WORLD) 
#define BCAST_BUF(buf, bufSize,type,who)  MPI_Bcast(buf, bufSize, type, who,MPI_COMM_WORLD )


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
#define VOLATILE_PAR volatile 
#define MASTER_P (tid == 0)
#define ASSIGN_INT(x,y) (x = y)
#define ASSIGN_BUF(x,y) (x = y)
#define ASSIGN_BUF_DBL(x,y) (x = y)
#define ASSIGN_DBL(x,y) (x = y)
#define ASSIGN_DBLS(tar,src,length) memmove(tar, src, length * sizeof(double))
#define DOUBLE double 	/* just rededining that to make the source code less confusing */
#define ASSIGN_GATHER(tar,src,length,type,tid) (memcpy((tar) + (tid) * (length) ,src, length * sizeof(type)))
#define SEND_BUF(buf, bufSize, type) 
#define RECV_BUF(buf, bufSize, type) 
#define BCAST_BUF(buf, bufSize,type,who)  
#endif


#endif	/* end include guard  */
