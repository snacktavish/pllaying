#ifndef MPI_HELPERS_H
#define MPI_HELPERS_H

int* addIntToBuf(int* buf, int *toAdd); 
int* popIntFromBuf(int *buf, int *result); 
void defineTraversalInfoMPI(MPI_Datatype *result); 

#endif
