#include "axml.h"

#ifdef _FINE_GRAIN_MPI
/* 
   This seems to be a very safe method to define your own mpi
   datatypes (often there are problems with padding). But it is not
   entirely for the weak of heart...
*/


int* addIntToBuf(int* buf, int *toAdd)
{
  *buf  = *toAdd; 
  return buf; 
}

int* popIntFromBuf(int *buf, int *result)
{
  *result = *buf; 
  return buf; 
}


/* :TODO: if we really want to overdo it, this could be templated => or use defines instead  */ 
double* addDblToBuf(double* buf, double *toAdd)
{
  *buf  = *toAdd; 
  return buf; 
}

double* popDblFromBuf(double *buf, double *result)
{
  *result = *buf; 
  return buf; 
}

#define ELEMS_IN_TRAV_INFO  9
void defineTraversalInfoMPI(MPI_Datatype *result)
{
/* :TODO: I guess, we have to free the mpi datatype later */
/* :TODO: in the best case, add all defined datatypes to an array in tree */

  int i ; 
  MPI_Aint base; 
  int blocklen[ELEMS_IN_TRAV_INFO+1] = {1, 1, 1, 1, NUM_BRANCHES, NUM_BRANCHES, 1,1,1,1}; 
  MPI_Aint disp[ELEMS_IN_TRAV_INFO+1];
  MPI_Datatype type[ELEMS_IN_TRAV_INFO+1] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_UB}; 
  traversalInfo desc[2]; 

  MPI_Address( desc, disp);
  MPI_Address( &(desc[0].pNumber), disp + 1 );
  MPI_Address( &(desc[0].qNumber), disp + 2 );  
  MPI_Address( &(desc[0].rNumber), disp + 3); 
  MPI_Address( desc[0].qz, disp + 4 );
  MPI_Address( desc[0].rz, disp + 5 );
  MPI_Address( &(desc[0].slot_p), disp + 6);
  MPI_Address( &(desc[0].slot_q), disp + 7);
  MPI_Address( &(desc[0].slot_r), disp + 8);
  MPI_Address( desc + 1, disp + 9);

  base = disp[0]; 
  for(i = 0; i < ELEMS_IN_TRAV_INFO+1; ++i)
    disp[i] -= base;

  MPI_Type_struct( ELEMS_IN_TRAV_INFO+1 , blocklen, disp, type, result);
  MPI_Type_commit(result);
}

#endif
