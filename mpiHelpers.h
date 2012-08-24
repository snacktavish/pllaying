#ifndef MPI_HELPERS_H
#define MPI_HELPERS_H

int* addIntToBuf(int* buf, int *toAdd); 
int* popIntFromBuf(int *buf, int *result); 
double* addDblToBuf(double* buf, double *toAdd); 
double* popDblFromBuf(double *buf, double *result); 
void defineTraversalInfoMPI(MPI_Datatype *result); 

#endif
