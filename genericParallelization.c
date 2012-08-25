#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include "genericParallelization.h"
#include "axml.h"


/* 
   GENERAL :TODO:
   unified system for tid/processId , numberOfThreads / processes
*/



/* 
 * GENERAL QUESTIONS: 
 * when to communicate rateCategories, wr, w2r 
 */

extern unsigned int* mask32; 
extern volatile int jobCycle; 
extern volatile int threadJob; 

void initializePartitionData(tree *localTree);
void initMemorySavingAndRecom(tree *tr); 

extern double *globalResult; 
extern volatile char *barrierBuffer;
void pinToCore(int tid);


/*****************/
/* MPI specific  */
/*****************/
#ifdef _FINE_GRAIN_MPI

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




/* 
   This seems to be a very safe method to define your own mpi
   datatypes (often there are problems with padding). But it is not
   entirely for the weak of heart...
*/

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


/********************/
/* PTHREAD-SPECIFIC */
/********************/
#ifdef _USE_PTHREADS

#ifndef _PORTABLE_PTHREADS
void pinToCore(int tid)
{
  cpu_set_t cpuset;

  CPU_ZERO(&cpuset);    
  CPU_SET(tid, &cpuset);

  if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
    {
      printBothOpen("\n\nThere was a problem finding a physical core for thread number %d to run on.\n", tid);
      printBothOpen("Probably this happend because you are trying to run more threads than you have cores available,\n");
      printBothOpen("which is a thing you should never ever do again, good bye .... \n\n");
      assert(0);
    }
}
#endif

void startPthreads(tree *tr)
{
  pthread_t *threads;
  pthread_attr_t attr;
  int rc, t;
  threadData *tData;

  jobCycle        = 0;
  threadJob       = 0;

  printf("\nThis is the RAxML Master Pthread\n");

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  threads    = (pthread_t *)malloc((size_t)tr->numberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc((size_t)tr->numberOfThreads * sizeof(threadData));

  barrierBuffer            = (volatile char *)  malloc(sizeof(volatile char)   *  (size_t)tr->numberOfThreads);

  for(t = 0; t < tr->numberOfThreads; t++)
    barrierBuffer[t] = 0;

  for(t = 1; t < tr->numberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
}
#endif



/* function that computes per-site log likelihoods in pthreads */

void perSiteLogLikelihoodsPthreads(tree *tr, double *lhs)
{
  int n = tr->numberOfThreads, 
    tid = tr->threadID; 
  size_t 
    model, 
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      size_t 
	localIndex = 0;

      /* decide if this partition is handled by the thread when -Q is ativated 
	 or when -Q is not activated figure out which sites have been assigned to the 
	 current thread */

      boolean 
	execute = ((tr->manyPartitions && isThisMyPartition(tr, tid, model)) || (!tr->manyPartitions));

      /* if the entire partition has been assigned to this thread (-Q) or if -Q is not activated 
	 we need to compute some per-site log likelihoods with thread tid for this partition */

      if(execute)
	for(i = tr->partitionData[model].lower;  i < tr->partitionData[model].upper; i++)
	  {
	    /* if -Q is active we compute all per-site log likelihoods for the partition,
	       othwerise we only compute those that have been assigned to thread tid 
	       using the cyclic distribution scheme */

	    if(tr->manyPartitions || (i % n == tid))
	      {
		double 
		  l;

		/* now compute the per-site log likelihood at the current site */

		switch(tr->rateHetModel)
		  {
		  case CAT:
		    l = evaluatePartialGeneric (tr, localIndex, tr->partitionData[model].perSiteRates[tr->partitionData[model].rateCategory[localIndex]], model);
		    break;
		  case GAMMA:
		    l = evaluatePartialGeneric (tr, localIndex, 1.0, model);
		    break;
		  default:
		    assert(0);
		  }

		/* store it in an array that is local in memory to the current thread,
		   see function collectDouble() in axml.c for understanding how we then collect these 
		   values stored in local arrays from the threads */

		lhs[i] = l;

		localIndex++;
	      }
	  }
    }
}



/* TODO: move to separate file?  */
boolean isThisMyPartition(tree *localTree, int tid, int model)
{ 
  if(localTree->partitionAssignment[model] == tid)
    return TRUE;
  else
    return FALSE;
}

void computeFractionMany(tree *localTree, int tid)
{
  int
    sites = 0;

  int   
    model;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      if(isThisMyPartition(localTree, tid, model))
	{	 
	  localTree->partitionData[model].width = localTree->partitionData[model].upper - localTree->partitionData[model].lower;
	  sites += localTree->partitionData[model].width;
	}
      else       	  
	localTree->partitionData[model].width = 0;       
    }


}



void computeFraction(tree *localTree, int tid, int n)
{
  int
    i,
    model;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      int width = 0;

      for(i = localTree->partitionData[model].lower; i < localTree->partitionData[model].upper; i++)
	if(i % n == tid)
	  width++;

      localTree->partitionData[model].width = width;
    }
}




static int partCompare(const void *p1, const void *p2)
{
  partitionType 
    *rc1 = (partitionType *)p1,
    *rc2 = (partitionType *)p2;

  int 
    i = rc1->partitionLength,
    j = rc2->partitionLength;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


/* 
   tr->manyPartitions is set to TRUE if the user has indicated via -Q that there are substantially more partitions 
   than threads/cores available. In that case we do not distribute sites from each partition in a cyclic fashion to the cores 
   , but distribute entire partitions to cores. 
   Achieving a good balance of alignment sites to cores boils down to the multi-processor scheduling problem known from theoretical comp. sci.
   which is NP-complete.
   We have implemented very simple "standard" heuristics for solving the multiprocessor scheduling problem that turn out to work very well
   and are cheap to compute. 
*/

void multiprocessorScheduling(tree *tr, int tid)
{
  int 
    s,
    model,
    modelStates[2] = {4, 20},
    numberOfPartitions[2] = {0 , 0},
      arrayLength = sizeof(modelStates) / sizeof(int);

      /* check that we have not addedd any new models for data types with a different number of states
	 and forgot to update modelStates */

      tr->partitionAssignment = (int *)malloc((size_t)tr->NumberOfModels * sizeof(int));

      for(model = 0; model < tr->NumberOfModels; model++)
	{        
	  boolean 
	    exists = FALSE;

	  for(s = 0; s < arrayLength; s++)
	    {
	      exists = exists || (tr->partitionData[model].states == modelStates[s]);
	      if(tr->partitionData[model].states == modelStates[s])
		numberOfPartitions[s] += 1;
	    }

	  assert(exists);
	}

      if(tid == 0)
	printBothOpen("\nMulti-processor partition data distribution enabled (\"-Q\" option)\n");

      for(s = 0; s < arrayLength; s++)
	{
	  if(numberOfPartitions[s] > 0)
	    {
	      size_t   
		checkSum = 0,
		sum = 0;

	      int    
		i,
		k,
#ifndef _FINE_GRAIN_MPI
		n = tr->numberOfThreads,
#else
		n = processes,
#endif
		p = numberOfPartitions[s],    
		*assignments = (int *)calloc((size_t)n, sizeof(int));  

	      partitionType 
		*pt = (partitionType *)malloc(sizeof(partitionType) * (size_t)p);



	      for(i = 0, k = 0; i < tr->NumberOfModels; i++)
		{
		  if(tr->partitionData[i].states == modelStates[s])
		    {
		      pt[k].partitionNumber = i;
		      pt[k].partitionLength = tr->partitionData[i].upper - tr->partitionData[i].lower;
		      sum += (size_t)pt[k].partitionLength;
		      k++;
		    }
		}

	      assert(k == p);

	      qsort(pt, p, sizeof(partitionType), partCompare);    

	      for(i = 0; i < p; i++)
		{
		  int 
		    k, 
		    min = INT_MAX,
		    minIndex = -1;

		  for(k = 0; k < n; k++)	
		    if(assignments[k] < min)
		      {
			min = assignments[k];
			minIndex = k;
		      }

		  assert(minIndex >= 0);

		  assignments[minIndex] +=  pt[i].partitionLength;
		  assert(pt[i].partitionNumber >= 0 && pt[i].partitionNumber < tr->NumberOfModels);
		  tr->partitionAssignment[pt[i].partitionNumber] = minIndex;
		}

	      if(tid == 0)
		{
		  for(i = 0; i < n; i++)	       
		    printBothOpen("Process %d has %d sites for %d state model \n", i, assignments[i], modelStates[s]);		  		

		  printBothOpen("\n");
		}

	      for(i = 0; i < n; i++)
		checkSum += (size_t)assignments[i];

	      assert(sum == checkSum);

	      free(assignments);
	      free(pt);
	    }
	}
}



/* 
   here again we collect the first and secdon derivatives from the various threads 
   and sum them up. It's similar to what we do in evaluateGeneric() 
   with the only difference that we have to collect two values (firsrt and second derivative) 
   instead of onyly one (the log likelihood

   Note: operates on global reduction buffers 
*/
void branchLength_parallelReduce(tree *tr, double *dlnLdlz,  double *d2lnLdlz2 ) 
{
  /* only the master executes this  */
  assert(tr->threadID == 0); 
  
  int b; 
  int t; 
  for(b = 0; b < tr->numBranches; ++b)
    {
      dlnLdlz[b] = 0; 
      d2lnLdlz2[b] = 0; 

      for(t = 0; t < tr->numberOfThreads; ++t)
	{
	  dlnLdlz[b] += globalResult[t * tr->numBranches * 2 + b ]; 
	  d2lnLdlz2[b] += globalResult[t * tr->numBranches * 2 + tr->numBranches + b]; 
	}
    } 
}



/* 
   reads from buffer or writes rates (doubles) into buffer 
   returns number  of elems written / read 
   countOnly: simply return the number of elements 
*/
static int doublesToBuffer(double *buf, double *srcTar, tree *tr, int n, int tid, boolean read, boolean countOnly)
{
  int 
    model,
    i;
  double 
    *initPtr = buf; 

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(tr, tid, model))	
	    for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i]; 
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++;
	      }	  
	}
      
      else
	{
	  for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)	    
	    if(i % n == tid)
	      {
		if(NOT countOnly)
		  {
		    if(read)
		      *buf = srcTar[i];
		    else 
		      srcTar[i] = *buf; 
		  }
		buf++; 
	      }
	}
    }
  
  return buf - initPtr; 
}


void gatherAndRedistributeDoubles(double *local, double *global, tree *tr, int n, int tid)
{  
  double 
    buf[tr->originalCrunchedLength]; 
  
  for(int i = 0; i < n; ++i )
    {      
#ifdef _USE_PTHREADS
      if(i != tid) 		/* for pthreads, only operate on private rates */
	continue;
#endif

      /* extract rates into buffer (for relevant thread) or just get the number of rates for this thread */
      int numberOfRates = doublesToBuffer(buf, local,tr, n , i, TRUE, i != tid); 
      
      /* each mpi worker broadcasts its rates (pthreads: no effect) */
      BCAST_BUF(buf, numberOfRates, MPI_DOUBLE, i); 

      /* now each process received the rates and writes it back into its structure */
      int assertNum = doublesToBuffer(buf, global, tr, n, i, FALSE, FALSE); 
      assert(assertNum == numberOfRates); 
    }  
}




static void collectDouble(double *dst, double *src, tree *tr, int n, int tid)
{
  double 
    resultBuf[tr->originalCrunchedLength],
    buf[tr->originalCrunchedLength]; 
  int
    assertNum, 
    displacements[tr->numberOfThreads]; 


  /* gather own values into buffer  */
  int numberCollected = doublesToBuffer(buf, src, tr,n,tid,TRUE, FALSE);   

#ifdef _FINE_GRAIN_MPI 
  /* this communicates all the values to the master */
  
  int numberPerWorker[tr->numberOfThreads];   
  if(MASTER_P)			/* master counts number to receive, receives and writes back */
    {
      for(int i = 0; i < n; ++i)
	{
	  numberPerWorker[i] = doublesToBuffer(buf,src,tr,n,i,FALSE, TRUE); 
	  displacements[i] = i == 0 ? 0 : displacements[i-1] + numberPerWorker[i-1]; 
	  printf("numPerW = %d\tdispl=%d\n", numberPerWorker[i], displacements[i]); 
	}
      
      MPI_Gatherv(buf, numberCollected, MPI_DOUBLE,
		  resultBuf, numberPerWorker, displacements,  MPI_DOUBLE,
		  0, MPI_COMM_WORLD); 

      printf("\n\nMASTER received"); 
      for(int i = 0; i < tr->originalCrunchedLength; ++i)
	printf("%f,", resultBuf[i]);
      printf("\n\n"); 

      double *bufPtr = resultBuf; 
      for(int i = 0 ; i < n; ++i)
	{
	  int numberWritten = doublesToBuffer(bufPtr, dst,tr,n,i, FALSE, FALSE); 
	  bufPtr += numberWritten; 
	  assertNum += numberWritten; 
	}    
      
      assert(assertNum == tr->originalCrunchedLength);
    }
  else 				/* workers only send their buffer   */
    MPI_Gatherv(buf, numberCollected, MPI_DOUBLE, resultBuf, numberPerWorker, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);   
#else 
  /* pthread only writes to global space  */  
  assertNum = doublesToBuffer(buf, dst,tr,n,tid, FALSE, FALSE);
  assert(assertNum == numberCollected); 
#endif
}


static void broadCastAlpha(tree *localTree, tree *tr, int tid)
{
  int 
    model; 

#ifdef _LOCAL_DISCRETIZATION 
  assert(0); 
  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
      makeGammaCats(localTree->partitionData[model].alpha, localTree->partitionData[model].gammaRates, 4, tr->useMedian); 
    }   
#else
  
  int
    i,
    bufSize = localTree->NumberOfModels * 4; 

  double bufDbl[bufSize], 
    *bufPtrDbl = bufDbl;   

  RECV_BUF(bufDbl, bufSize, MPI_DOUBLE); 

  for(model = 0; model < localTree->NumberOfModels; model++)
    for(i = 0; i < 4; ++i)
      ASSIGN_BUF_DBL(localTree->partitionData[model].gammaRates[i], tr->partitionData[model].gammaRates[i]); 
  
  SEND_BUF(bufDbl, bufSize, MPI_DOUBLE);  
#endif 
}


static void broadCastRates(tree *localTree, tree *tr, int tid)
{
  int 
    model;
#ifdef _LOCAL_DISCRETIZATION 
  assert(0); 
  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
      
      if(tid > 0)
	memcpy(localTree->partitionData[model].substRates,        tr->partitionData[model].substRates, pl->substRatesLength * sizeof(double));
      
      initReversibleGTR(localTree, model);     
    } 
#else

  /* determine size of buffer needed first */
  int bufSize = 0; 
#ifdef _FINE_GRAIN_MPI

  for(model = 0; model < localTree->NumberOfModels; ++model )
    {	  
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model])); /* this is constant, isnt it?  */
      bufSize += pl->eignLength + pl->evLength + pl->eiLength + pl->tipVectorLength; 
    }
#endif      
  
  /* 
     the usual game: MPI version writes everything into a buffer that
     is broadcasted; pthreads version uses direct assignment. 

     TODO memcpy still is more efficient
  */
  
  double bufDbl[bufSize],
    *bufPtrDbl  = bufDbl;
  RECV_BUF(bufDbl, bufSize, MPI_DOUBLE);
  int i ; 

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model])); /* this is constant, isnt it?  */

      for(i = 0; i < pl->eignLength; ++i)
	ASSIGN_BUF_DBL(localTree->partitionData[model].EIGN[i], tr->partitionData[model].EIGN[i]); 
      for(i = 0; i < pl->evLength; ++i)
	ASSIGN_BUF_DBL(localTree->partitionData[model].EV[i],tr->partitionData[model].EV[i]);
      for(i = 0; i  < pl->eiLength; ++i)
	ASSIGN_BUF_DBL(localTree->partitionData[model].EI[i], tr->partitionData[model].EI[i]);
      for(i = 0; i < pl->tipVectorLength; ++i)
	ASSIGN_BUF_DBL(localTree->partitionData[model].tipVector[i],   tr->partitionData[model].tipVector[i]);
    }
  SEND_BUF(bufDbl, bufSize, MPI_DOUBLE); 
#endif
}


static void reduceEvaluateIterative(tree *localTree, int tid)
{
  int model;

  evaluateIterative(localTree);

  /* when this is done we need to write the per-thread log likelihood to the 
     global reduction buffer. Tid is the thread ID, hence thread 0 will write its 
     results to reductionBuffer[0] thread 1 to reductionBuffer[1] etc.

     the actual sum over the entries in the reduction buffer will then be computed 
     by the master thread which ensures that the sum is determinsitic */


  /* 
     aberer: i implemented this as a mpi_gather operation into this buffer, 
     pthreads version emulates this gather; 
     master takes care of the reduction; 
  */

  double buf[localTree->NumberOfModels]; 
  for(model = 0; model < localTree->NumberOfModels; ++model)
    buf[model] = localTree->perPartitionLH[model];        

  ASSIGN_GATHER(globalResult, buf, localTree->NumberOfModels, DOUBLE, tid); 
}


/* the one below is a hack we are re-assigning the local pointer to the global one
   the memcpy version below is just for testing and preparing the
   fine-grained MPI BlueGene version */
/* TODO: we should reset this at some point, the excplicit copy is just done for testing */
inline static void broadcastTraversalInfo(tree *localTree, tree *tr)
{

  ASSIGN_INT(localTree->td[0].functionType,           tr->td[0].functionType);
  ASSIGN_INT( localTree->td[0].count,                 tr->td[0].count) ;
  ASSIGN_INT(localTree->td[0].traversalHasChanged,    tr->td[0].traversalHasChanged);

#ifdef _USE_PTHREADS
  /* memcpy -> memmove (see ticket #43). This function is sometimes called with localTree == tr,
   * in which case some memcpy implementations can corrupt the buffers.
   */ 

  memmove(localTree->td[0].executeModel,    tr->td[0].executeModel,    sizeof(boolean) * localTree->NumberOfModels);
  memmove(localTree->td[0].parameterValues, tr->td[0].parameterValues, sizeof(double) * localTree->NumberOfModels);
  
  if(localTree->td[0].traversalHasChanged)
    memmove(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));

#else  /* MPI */
  localTree = tr; 

  if(localTree->td[0].functionType != THREAD_INIT_PARTITION)				/* :KLUDGE: is not initialized in this case   */
    {
      MPI_Bcast(localTree->td[0].executeModel, localTree->NumberOfModels, MPI_INT, 0,MPI_COMM_WORLD); 
      MPI_Bcast(localTree->td[0].parameterValues, localTree->NumberOfModels, MPI_DOUBLE, 0,MPI_COMM_WORLD); 
      if(localTree->td[0].traversalHasChanged)
	{
	  /* define the datatype and broadcast */
	  MPI_Datatype trav_MPI; 
	  defineTraversalInfoMPI(&trav_MPI); 
	  MPI_Bcast(localTree->td[0].ti, localTree->td[0].count, trav_MPI, 0, MPI_COMM_WORLD ); 
	}
    }
#endif
}




void printParallelDebugInfo(int type, int tid )
{
  switch(type)  
    {
    case  THREAD_NEWVIEW: 
      printf("[%d] working on  THREAD_NEWVIEW\n ", tid); 
      break; 
    case THREAD_EVALUATE: 
      printf("[%d] working on  THREAD_EVALUATE\n ", tid); 
      break; 
    case THREAD_MAKENEWZ: 
      printf("[%d] working on  THREAD_MAKENEWZ\n ", tid); 
      break; 
    case THREAD_MAKENEWZ_FIRST: 
      printf("[%d] working on  THREAD_MAKENEWZ_FIRST\n ", tid); 
      break; 
    case THREAD_RATE_CATS: 
      printf("[%d] working on  THREAD_RATE_CATS\n ", tid); 
      break; 
    case THREAD_COPY_RATE_CATS: 
      printf("[%d] working on  THREAD_COPY_RATE_CATS\n ", tid); 
      break; 
    case THREAD_COPY_INIT_MODEL: 
      printf("[%d] working on  THREAD_COPY_INIT_MODEL\n ", tid); 
      break; 
    case THREAD_INIT_PARTITION: 
      printf("[%d] working on  THREAD_INIT_PARTITION\n ", tid); 
      break; 
    case THREAD_OPT_ALPHA: 
      printf("[%d] working on  THREAD_OPT_ALPHA\n ", tid); 
      break; 
    case THREAD_OPT_RATE: 
      printf("[%d] working on  THREAD_OPT_RATE\n ", tid); 
      break; 
    case THREAD_COPY_ALPHA: 
      printf("[%d] working on  THREAD_COPY_ALPHA\n ", tid); 
      break; 
    case THREAD_COPY_RATES: 
      printf("[%d] working on  THREAD_COPY_RATES\n ", tid); 
      break; 
    case THREAD_PER_SITE_LIKELIHOODS: 
      printf("[%d] working on  THREAD_PER_SITE_LIKELIHOODS\n ", tid); 
      break; 
    case THREAD_NEWVIEW_ANCESTRAL: 
      printf("[%d] working on  THREAD_NEWVIEW_ANCESTRAL\n ", tid); 
      break; 
    case THREAD_GATHER_ANCESTRAL: 
      printf("[%d] working on  THREAD_GATHER_ANCESTRAL\n ", tid); 
      break; 
    case THREAD_EXIT_GRACEFULLY: 
      printf("[%d] working on  THREAD_EXIT_GRACEFULLY\n ", tid); 
      break; 
    default: assert(0); 
    }  
}


/* this function here handles all parallel regions in the Pthreads version, when we enter 
   this function masterBarrier() has ben called by the master thread from within the sequential 
   part of the program, tr is the tree at the master thread, localTree the tree at the worker threads

   While this is not necessary, adress spaces of threads are indeed separated for easier transition to 
   a distributed memory paradigm 
*/

boolean execFunction(tree *tr, tree *localTree, int tid, int n)
{
  int
    i,
    model,
    localCounter;

#ifdef _USE_PTHREADS
  /* some stuff associated with the barrier implementation using Pthreads and busy wait */
  int currentJob = threadJob >> 16;
#endif

  /* here the master sends and all threads/processes receive the traversal descriptor */
  broadcastTraversalInfo(localTree, tr);
  
#ifdef _USE_PTHREADS
  /* make sure that nothing is going wrong */
  assert(currentJob == localTree->td[0].functionType);
#else   
  localTree = tr; 
  int currentJob = localTree->td[0].functionType; 
#endif
#ifdef DEBUG_PARALLEL
  /* printf("process %d: current job is %d\n", tid, currentJob);  */
  printParallelDebugInfo(currentJob, tid);
#endif  

  switch(currentJob)
    { 
    case THREAD_NEWVIEW:      
      /* just a newview on the fraction of sites that have been assigned to this thread */

      newviewIterative(localTree, 0);
      break;     
    case THREAD_EVALUATE: 
      reduceEvaluateIterative(localTree, tid); 
      break;	
    case THREAD_MAKENEWZ_FIRST:

      /* this is the first call from within makenewz that requires getting the likelihood vectors to the left and 
         right of the branch via newview and doing some precomputations.
	 
         For details see comments in makenewzGenericSpecial.c 
      */
    case  THREAD_MAKENEWZ:
      {	
	double
	  dlnLdlz[NUM_BRANCHES],
	  d2lnLdlz2[NUM_BRANCHES]; 

	if(localTree->td[0].functionType == THREAD_MAKENEWZ_FIRST)
	  makenewzIterative(localTree);	
	execCore(localTree, dlnLdlz, d2lnLdlz2);

	/* gather the first and second derivatives that have been written by each thread */
	/* as for evaluate above, the final sum over the derivatives will be computed by the 
	   master thread in its sequential part of the code */

	/* MPI: implemented as a gather again, pthreads: just buffer copying */	
	double buf[ 2 * localTree->numBranches]; 
	memcpy( buf, dlnLdlz, localTree->numBranches * sizeof(double) ); 
	memcpy(buf + localTree->numBranches, d2lnLdlz2, localTree->numBranches * sizeof(double)); 
	
	ASSIGN_GATHER(globalResult, buf,  2 * localTree->numBranches, DOUBLE, tid); 

      }

      break;

    case THREAD_INIT_PARTITION:       
      
      /* broadcast data and initialize and allocate arrays in partitions */
      
      initializePartitions(tr, localTree, tid, n);
      
      break;          
    case THREAD_COPY_ALPHA: 
    case THREAD_OPT_ALPHA:
      /* this is when we have changed the alpha parameter, inducing a change in the discrete gamma rate categories.
	 this is called when we are optimizing or sampling (in the Bayesioan case) alpha parameter values */
      
      /* distribute the new discrete gamma rates to the threads */
      broadCastAlpha(localTree,tr, tid);

      /* compute the likelihood, note that this is always a full tree traversal ! */
      if(localTree->td[0].functionType == THREAD_OPT_ALPHA)
	reduceEvaluateIterative(localTree, tid);      

      break;           
    case THREAD_OPT_RATE:
    case THREAD_COPY_RATES:

      /* if we are optimizing the rates in the transition matrix Q this induces recomputing the eigenvector eigenvalue 
	 decomposition and the tipVector as well because of the special numerics in RAxML, the matrix of eigenvectors 
	 is "rotated" into the tip lookup table.

	 Hence if the sequantial part of the program that steers the Q matrix rate optimization has changed a rate we 
	 need to broadcast all eigenvectors, eigenvalues etc to each thread 
      */

      broadCastRates(localTree, tr, tid);     

      /* now evaluate the likelihood of the new Q matrix, this always requires a full tree traversal because the changes need
	 to be propagated throughout the entire tree */

      if(localTree->td[0].functionType == THREAD_OPT_RATE)
	reduceEvaluateIterative(localTree, tid); 

      break;                       
    case THREAD_COPY_INIT_MODEL:
      {
	/* 
	   need to copy base freqs before local per-thread/per-process Q matrix exponentiation 
	*/

#ifdef _LOCAL_DISCRETIZATION   
	assert(0); /* andre: did not implement this part, since it is not sure, if it persists   */

	if(tid > 0)
	  for(model = 0; model < localTree->NumberOfModels; model++)	    	     
	    {
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
	 
	      if(tid > 0)
		memcpy(localTree->partitionData[model].frequencies,        tr->partitionData[model].frequencies,        pl->frequenciesLength * sizeof(double));
	    }
#endif
	/* need to be very careful here ! THREAD_COPY_INIT_MODEL is also used when the program is restarted 
	   it is hence not sufficient to just initialize everything by the default values ! */

	broadCastRates(localTree, tr, tid); 
	broadCastAlpha(localTree, tr, tid); /* isnt that only executed when we are on gamma?  */

	/*
	  copy initial model parameters, the Q matrix and alpha are initially, when we start our likelihood search 
	  set to default values. 
	  Hence we need to copy all those values that are required for computing the likelihood 
	  with newview(), evaluate() and makenez() to the private memory of the threads 
	*/


	if( localTree->rateHetModel == CAT) /* TRICKY originally this should only be executed by workers  */
	  { 
	    int bufSize = 2 * localTree->originalCrunchedLength; 
	    double 
	      bufDbl[bufSize],
	      *bufPtrDbl = bufDbl; 	 

	    RECV_BUF(bufDbl, bufSize,MPI_DOUBLE); 

	    /* this should be local  */
	    for(model = 0; model < localTree->NumberOfModels; model++) 
	      localTree->partitionData[model].numberOfCategories      = tr->partitionData[model].numberOfCategories;	    


	    /* this is only relevant for the PSR model, we can worry about this later */
	    for(i = 0; i < localTree->originalCrunchedLength; ++i)
	      {
		ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]);
		ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	      }

	    SEND_BUF(bufDbl, bufSize, MPI_DOUBLE); 
	  }
      } 
      break;    
    case THREAD_RATE_CATS: 
      {
	/* this is for optimizing per-site rate categories under PSR, let's worry about this later */

	ASSIGN_DBL( localTree->lower_spacing,  tr->lower_spacing);
	ASSIGN_DBL( localTree->upper_spacing,  tr->upper_spacing);

	optRateCatPthreads(localTree, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

	gatherAndRedistributeDoubles(localTree->patrat, tr->patrat, localTree, n, tid); 
	gatherAndRedistributeDoubles(localTree->patratStored, tr->patratStored, localTree, n, tid); 
	gatherAndRedistributeDoubles(localTree->lhs, tr->lhs, localTree, n, tid); 

      }
      break;
    case THREAD_COPY_RATE_CATS:
      {
	/* 
	   this is invoked when we have changed the per-site rate category assignment
	   In essence it distributes the new per site rates to all threads 

	   The pthread-version here simply assigns everything as ought to
	   be. The MPI-version is configured to write to a buffer instead
	   and SEND (master) or RECV (workers) it.

	*/

	/* 
	   start of communication part 
	*/

	int i, 
	  buf[localTree->NumberOfModels],
	  assertCtr = 0, 
	  *bufPtr = buf, 
	  dblBufSize = 0; 
     
	RECV_BUF(buf, localTree->NumberOfModels, MPI_INT);

	for( model = 0; model < localTree->NumberOfModels; ++model)
	  {
	    ASSIGN_BUF(localTree->partitionData[model].numberOfCategories, tr->partitionData[model].numberOfCategories); 
	    dblBufSize += localTree->partitionData[model].numberOfCategories; 
	  }

	SEND_BUF(buf, localTree->NumberOfModels, MPI_INT); 


	dblBufSize += 2 * localTree->originalCrunchedLength; 
	double
	  bufDbl[dblBufSize], 
	  *bufPtrDbl = bufDbl;      

	RECV_BUF(bufDbl, dblBufSize, MPI_DOUBLE);      

	for(i = 0; i < localTree->originalCrunchedLength; ++i)
	  {	 
	    ASSIGN_BUF_DBL(localTree->patrat[i], tr->patrat[i]); 
	    ASSIGN_BUF_DBL(localTree->patratStored[i], tr->patratStored[i]); 
	  }

	for( model = 0; model < localTree->NumberOfModels; ++model)
	  for(i = 0; i < localTree->partitionData[model].numberOfCategories; i++)
	    ASSIGN_BUF_DBL(localTree->partitionData[model].perSiteRates[i], tr->partitionData[model].perSiteRates[i]);

	SEND_BUF(bufDbl, dblBufSize, MPI_DOUBLE); 




	/* lets test, if it is a good idea to send around the basic categories  */
#ifdef _FINE_GRAIN_MPI
	/* TODO: this is inefficient, need to clarify how much this is necessary   */
	MPI_Bcast(tr->rateCategory, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD); 
	MPI_Bcast(tr->wr2, tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
	MPI_Bcast(tr->wr, tr->originalCrunchedLength, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
#endif


	/* 
	   now re-assign values 
	*/
	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    if(localTree->manyPartitions)
	      {
		if(isThisMyPartition(localTree, tid, model))
		  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++, localCounter++)
		    {	     
		      localTree->partitionData[model].rateCategory[localCounter] = tr->rateCategory[i];
		      localTree->partitionData[model].wr[localCounter]             = tr->wr[i];
		      localTree->partitionData[model].wr2[localCounter]            = tr->wr2[i]; 
		    } 
	      }
	    else	  
	      {
		for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
		  {
		    if(i % n == tid)
		      {		 
			localTree->partitionData[model].rateCategory[localCounter] = tr->rateCategory[i]; 
			localTree->partitionData[model].wr[localCounter]             = tr->wr[i];
			localTree->partitionData[model].wr2[localCounter]            = tr->wr2[i];		 

			localCounter++;
		      }
		  }
	      }
	  }
      }
      break;
    case THREAD_PER_SITE_LIKELIHOODS:      
      {
	/* compute per-site log likelihoods for the sites/partitions 
	   that are handled by this thread */
	perSiteLogLikelihoodsPthreads(localTree, localTree->lhs);

	/* do a parallel gather operation, the threads will write their results 
	   into the global buffer tr->lhs that will then contain all per-site log likelihoods
	   in the proper order 
	*/

	collectDouble(tr->lhs,                localTree->lhs,                  localTree, n, tid);

	for(int i = 0; i < tr->originalCrunchedLength; ++i)
	  printf("%f,", tr->lhs[i]); 
	printf("\n"); 
	
      }
      break;
      /* check for errors */
    case THREAD_NEWVIEW_ANCESTRAL:       
      assert(0);
      break; 
    case THREAD_GATHER_ANCESTRAL:
      assert(0); 
      break; 
    case THREAD_EXIT_GRACEFULLY: 
#ifdef _FINE_GRAIN_MPI
      MPI_Finalize(); 
#endif
      return FALSE; 
      break; 
    default:
      printf("Job %d\n", currentJob);
      assert(0);
    }

  return TRUE; 
}


void *likelihoodThread(void *tData)
{
  threadData *td = (threadData*)tData;
  tree
    *tr = td->tr;

#ifdef _USE_PTHREADS
  tree *localTree = calloc(1,sizeof(tree )); 
  
  int
    myCycle = 0;

  const int 
    n = td->tr->numberOfThreads,
    tid = td->threadNumber;

#ifndef _PORTABLE_PTHREADS
  pinToCore(tid);
#endif

  printf("\nThis is RAxML Worker Pthread Number: %d\n", tid);

  while(1)
    {

      while (myCycle == threadJob);
      myCycle = threadJob;

      execFunction(tr, localTree, tid, n);

      barrierBuffer[tid] = 1;     
    }
#else 
  const int
    n = processes, 
    tid = ((threadData*)tData)->threadNumber;

  printf("\nThis is RAxML Worker Process Number: %d\n", tid);

  while(execFunction(tr,tr, tid,n)); 
#endif

  return (void*)NULL;
}





/* 
   This is master specific code called once the barrier is
   passed. Stuff such as reduction operations.  If we execute this
   here, we can keep the code mostly free from parallel -specific
   code.
*/
void masterPostBarrier(int jobType, tree *tr)
{
  int 
    i, 
    model; 

  assert(tr->threadID == 0); 
  
  switch(jobType)
    {
    case THREAD_EVALUATE: 
      {      
	for(model = 0; model < tr->NumberOfModels; model++)
	  { 
	    volatile double 
	      partitionResult = 0.0;  
	    
	    for(i = 0, partitionResult = 0.0; i < tr->numberOfThreads; i++)          	      
	      partitionResult += globalResult[i * tr->NumberOfModels + model]; 
	    
	    tr->perPartitionLH[model] = partitionResult;
#ifdef DEBUG_PARALLEL
	    printf("[%d] result for partition %d => %f\n", model, partitionResult); 
#endif
	  }
      }
    case THREAD_OPT_ALPHA: 
      {
	int j;
	volatile double partitionResult;	
	volatile double result = 0.0;

	for(j = 0; j < tr->NumberOfModels; j++)
	  {
	    for(i = 0, partitionResult = 0.0; i < tr->numberOfThreads; i++) 
	      partitionResult += globalResult[i * tr->NumberOfModels + j];
	    
	    result +=  partitionResult;
	    tr->perPartitionLH[j] = partitionResult;
	  }
	break; 
      }      
    case THREAD_OPT_RATE: 
      {
	volatile double result;	

	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < tr->numberOfThreads; i++)    	  
	      result += globalResult[i]; 		

	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    int j;
	    volatile double partitionResult;
	
	    result = 0.0;

	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < tr->numberOfThreads; i++)          	      
		  partitionResult += globalResult[i * tr->NumberOfModels + j]; 

		result +=  partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }	
	break; 
      }
    case THREAD_PER_SITE_LIKELIHOODS:
      {
	/* now just compute the sum over per-site log likelihoods for error checking */      
	double accumulatedPerSiteLikelihood = 0.; 
	for(i = 0; i < tr->originalCrunchedLength; i++)
	  accumulatedPerSiteLikelihood += tr->lhs[i];

	printf("RESULT: %f\t%f", tr->likelihood, accumulatedPerSiteLikelihood); 
	assert(ABS(tr->likelihood - accumulatedPerSiteLikelihood) < 0.00001);
      }
      break;
    } 
}


void masterBarrier(int jobType, tree *tr)
{

#ifdef _USE_PTHREADS
  const int 
    n = tr->numberOfThreads;

  tr->td[0].functionType = jobType;

  jobCycle = !jobCycle;
  threadJob = (jobType << 16) + jobCycle;

  execFunction(tr, tr, 0, n);
  
  int 
    i, 
    sum;

  do
    {
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);  

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
#else 
  tr->td[0].functionType = jobType; 
  execFunction(tr,tr,0,processes);
#endif

  /* code executed by the master, once the barrier is crossed */
  masterPostBarrier(jobType, tr);
}


#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

/* encapsulated this, s.t. it becomes more clear, that the pthread-master must not execute this */
static void assignAndInitPart1(tree *localTree, tree *tr, int *tid)
{
  size_t
    model; 
  int
    totalLength = 0; 

#ifdef _USE_PTHREADS

  localTree->threadID = *tid; 
  printf("my id is %d\n", *tid); 
  assert(localTree != tr);
  localTree->numberOfThreads = tr->numberOfThreads;
#else  /* => MPI */
  int 
    assertCtr = 0;

  *tid = processID; 
  localTree->threadID = processID; 
  tr->numberOfThreads = processes;
  int bufSize = 8 + tr->NumberOfModels * 8;
#ifdef _LOCAL_DISCRETIZATION
  bufSize++;
#endif
  int buf[bufSize]; 
  int *bufPtr = buf; 
  if(NOT MASTER_P)		/* :TRICKY: workers receive buffer here */
    MPI_Bcast(buf, bufSize, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  ASSIGN_BUF( localTree->useRecom,                  tr->useRecom);
  ASSIGN_BUF( localTree->rateHetModel,              tr->rateHetModel);
  ASSIGN_BUF( localTree->useMedian,                 tr->useMedian); 
  ASSIGN_BUF( localTree->saveMemory,                tr->saveMemory);
  ASSIGN_BUF( localTree->maxCategories,             tr->maxCategories);
  ASSIGN_BUF( localTree->originalCrunchedLength,    tr->originalCrunchedLength);
  ASSIGN_BUF( localTree->mxtips,                    tr->mxtips);
  ASSIGN_BUF( localTree->numBranches,               tr->numBranches);
#ifdef _LOCAL_DISCRETIZATION
  ASSIGN_BUF( localTree->aliaswgt,                  tr->aliaswgt);
#endif    

  localTree->td[0].count = 0; 

  if(NOT MASTER_P)
    {
      localTree->lhs                     = (double*)malloc(sizeof(double)   * (size_t)localTree->originalCrunchedLength);     
      localTree->perPartitionLH          = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);     
      localTree->fracchanges             = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);
      localTree->partitionContributions  = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);
      localTree->partitionData           = (pInfo*)malloc(sizeof(pInfo)     * (size_t)localTree->NumberOfModels);
      localTree->td[0].ti              = (traversalInfo *)malloc(sizeof(traversalInfo) * (size_t)localTree->mxtips);
      localTree->td[0].executeModel    = (boolean *)malloc(sizeof(boolean) * (size_t)localTree->NumberOfModels);
      localTree->td[0].parameterValues = (double *)malloc(sizeof(double) * (size_t)localTree->NumberOfModels);
      localTree->patrat       = (double*)malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);
      localTree->patratStored = (double*)malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);            
    }
  
  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      ASSIGN_BUF(      localTree->partitionData[model].numberOfCategories,     tr->partitionData[model].numberOfCategories);
      ASSIGN_BUF(      localTree->partitionData[model].states,                 tr->partitionData[model].states);
      ASSIGN_BUF(      localTree->partitionData[model].maxTipStates ,          tr->partitionData[model].maxTipStates);
      ASSIGN_BUF(      localTree->partitionData[model].dataType ,              tr->partitionData[model].dataType);
      ASSIGN_BUF(      localTree->partitionData[model].protModels ,            tr->partitionData[model].protModels);
      ASSIGN_BUF(      localTree->partitionData[model].protFreqs ,             tr->partitionData[model].protFreqs);
      ASSIGN_BUF(      localTree->partitionData[model].lower ,                 tr->partitionData[model].lower);
      ASSIGN_BUF(      localTree->partitionData[model].upper ,                 tr->partitionData[model].upper); 

      localTree->perPartitionLH[model]                      = 0.0;
      totalLength += (localTree->partitionData[model].upper -  localTree->partitionData[model].lower);
    }

#ifdef _FINE_GRAIN_MPI
  assert(assertCtr == bufSize); /* will fail, if anybody changes structures, s.t. more variables need to be sent around */
  if(MASTER_P)			/* :TRICKY: master broadcasts buffer here  */
    MPI_Bcast(buf, bufSize, MPI_INT, 0, MPI_COMM_WORLD);
#endif

  assert(totalLength == localTree->originalCrunchedLength);

  ASSIGN_DBL(localTree->vectorRecomFraction, tr->vectorRecomFraction); 
}
#endif


void distributeYVectors(tree *localTree, tree *tr)
{
  size_t 
    i,
    n = localTree->numberOfThreads,
    globalCounter = 0,
    localCounter = 0,
    model = 0, 
    j; 
  int tid = localTree->threadID; 
  

  /* distribute the y-vectors */
  for(j = 1 ; j <= (size_t)localTree->mxtips; j++)	
    {
#ifdef _FINE_GRAIN_MPI
      unsigned char yBuf[tr->originalCrunchedLength]; 	  
      if(MASTER_P)
	memcpy(yBuf, tr->yVector[j], tr->originalCrunchedLength * sizeof(unsigned char));
      MPI_Bcast(  yBuf, tr->originalCrunchedLength, MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD); 
#endif	  

      for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
	{
	  if(tr->manyPartitions)
	    {
	      if(isThisMyPartition(localTree, tid, model))
		{
		  assert(localTree->partitionData[model].upper - localTree->partitionData[model].lower == localTree->partitionData[model].width);
		  for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++, localCounter++, globalCounter++)
#ifdef _USE_PTHREADS
		    localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		  localTree->partitionData[model].yVector[j][localCounter] = yBuf[globalCounter];
#endif


		}
	      else
		globalCounter += (localTree->partitionData[model].upper - localTree->partitionData[model].lower);
	    }
	  else 
	    {
	      for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++, globalCounter++)
		{
		  if(i % (size_t)n == (size_t)tid)
		    {
#ifdef _USE_PTHREADS
		      localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter];
#else 
		      localTree->partitionData[model].yVector[j][localCounter] = yBuf[globalCounter];
#endif
		      ++localCounter; 
		    }
		}	   
	    }
	}
    }
}




void distributeWeights(tree *localTree, tree *tr)
{
  int tid = localTree->threadID; 
  int n = localTree->numberOfThreads; 

  size_t     
    globalCounter = 0,
    i,
    localCounter  = 0,
    model; 



  /* distribute the weights  */
#ifdef _FINE_GRAIN_MPI 		/* need to broadcast a few things first */
  if(NOT MASTER_P)
    tr->aliaswgt = malloc(sizeof(int) * tr->originalCrunchedLength); 
  MPI_Bcast(tr->aliaswgt, tr->originalCrunchedLength, MPI_INT, 0, MPI_COMM_WORLD);      
#endif
  for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
    { 
      if(tr->manyPartitions)
	{
	  if(isThisMyPartition(localTree, tid, model))
	    {
	      assert(localTree->partitionData[model].upper - localTree->partitionData[model].lower == localTree->partitionData[model].width);
	      for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++, localCounter++, globalCounter++)
		localTree->partitionData[model].wgt[localCounter]          = tr->aliaswgt[globalCounter]; 
	    }
	  else
	    globalCounter += (localTree->partitionData[model].upper - localTree->partitionData[model].lower);
	}
      else 
	{ 
	  for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++, globalCounter++)
	    {
	      if(i % (size_t)n == (size_t)tid)
		localTree->partitionData[model].wgt[localCounter++]       = tr->aliaswgt[globalCounter]; 
	    }	   
	}
    }
}


void initializePartitionsMaster(tree *tr, tree *localTree, int tid, int n)
{ 
  size_t
    model;

#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))

  ASSIGN_INT(localTree->manyPartitions, tr->manyPartitions);
  ASSIGN_INT(localTree->NumberOfModels, tr->NumberOfModels); 

#ifdef _USE_PTHREADS
  if(MASTER_P)
#endif
    globalResult = calloc((size_t) tr->numberOfThreads * (size_t)tr->NumberOfModels * 2 ,sizeof(double)); 
#ifdef _USE_PTHREADS    
  else 
#endif
    assignAndInitPart1(localTree, tr , &tid); 
  
  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    localTree->partitionData[model].width        = 0;
  
  if(tr->manyPartitions)    
    {
      multiprocessorScheduling(localTree, tid);        
      computeFractionMany(localTree, tid);
    }
  else
    computeFraction(localTree, tid, n);

#else  /* sequential */
  assert(tr == localTree);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    assert(tr->partitionData[model].width == tr->partitionData[model].upper - tr->partitionData[model].lower);
#endif	   

  initializePartitionData(localTree); 

#if NOT (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  /* figure in tip sequence data per-site pattern weights */ 
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t
	j,
	lower = tr->partitionData[model].lower,
	width = tr->partitionData[model].upper - lower;  

      for(j = 1; j <= (size_t)tr->mxtips; j++)
	tr->partitionData[model].yVector[j] = &(tr->yVector[j][tr->partitionData[model].lower]);

      memcpy((void*)(&(tr->partitionData[model].wgt[0])),         (void*)(&(tr->aliaswgt[lower])),      sizeof(int) * width);            
    }  
#else
  {
    size_t 
      model,  
      j,
      i,
      globalCounter = 0,
      localCounter  = 0,
      countOffset,
      myLength = 0;

    for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
      myLength += localTree->partitionData[model].width;         

    /* assign local memory for storing sequence data */
    
    localTree->y_ptr = (unsigned char *)malloc(myLength * (size_t)(localTree->mxtips) * sizeof(unsigned char));
    assert(localTree->y_ptr != NULL);

    for(i = 0; i < (size_t)localTree->mxtips; i++)
      {
	for(model = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
	  {	    
	    localTree->partitionData[model].yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
	    countOffset +=  localTree->partitionData[model].width;
	  }
	assert(countOffset == myLength);
      }

    /* figure in data */

    distributeWeights(localTree, tr); 

    distributeYVectors(localTree, tr); 
  }
#endif

  initMemorySavingAndRecom(localTree);
}


/* interface to outside  */
void initializePartitions(tree *tr, tree *localTree, int tid, int n)
{
  initializePartitionsMaster(tr,localTree,tid,n); 
}



void initMemorySavingAndRecom(tree *tr)
{
#ifdef DEBUG_PARALLEL
  printf("process %d: entering memorySavingRecom()\n", tr->threadID); 
#endif

  tree
    *localTree = tr; 
  size_t model; 

  /* initialize gap bit vectors at tips when memory saving option is enabled */

  if(localTree->saveMemory)
    {
      for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
	{
	  int        
	    undetermined = getUndetermined(localTree->partitionData[model].dataType);

	  size_t
	    i,
	    j,
	    width =  localTree->partitionData[model].width;

	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(localTree->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(localTree->partitionData[model].yVector[j][i] == undetermined)
		    localTree->partitionData[model].gapVector[localTree->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }
  /* recom */
  if(localTree->useRecom)
    allocRecompVectorsInfo(localTree);
  else
    localTree->rvec = (recompVectors*)NULL;
  /* E recom */

#ifdef DEBUG_PARALLEL
  printf("process %d: leaving memorySavingRecom()\n", tr->threadID); 
#endif
}


/* mostly malloc calls and initialization  */
void initializePartitionData(tree *localTree)
{
  /* in ancestralVectorWidth we store the total length in bytes (!) of 
     one conditional likelihood array !
     we need to know this length such that in the pthreads version the master thread can actually 
     gather the scattered ancestral probabilities from the threads such that they can be printed to screen!
  */

  size_t 
    maxCategories = (size_t)localTree->maxCategories;

  size_t 
    ancestralVectorWidth = 0,
    model; 
  int tid  = localTree->threadID; 


  /* aberer: added this, we need it for mpi    */
  if(NOT MASTER_P)
    {
      localTree->rateCategory    = (int *)    calloc((size_t)localTree->originalCrunchedLength, sizeof(int));	    
      localTree->wr              = (double *) calloc((size_t)localTree->originalCrunchedLength, sizeof(double)); 
      localTree->wr2             = (double *) calloc((size_t)localTree->originalCrunchedLength, sizeof(double));   
    }
  /* end aberer */

  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      size_t 
	j,       
	width = localTree->partitionData[model].width;

      const partitionLengths 
	*pl = getPartitionLengths(&(localTree->partitionData[model]));

      localTree->partitionData[model].wr = (double *)malloc(sizeof(double) * width);
      localTree->partitionData[model].wr2 = (double *)malloc(sizeof(double) * width);     


      /* 
	 globalScaler needs to be 2 * localTree->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
      */

      localTree->partitionData[model].globalScaler    = (unsigned int *)calloc(2 *(size_t)localTree->mxtips, sizeof(unsigned int));  	         

      localTree->partitionData[model].left              = (double *)malloc_aligned((size_t)pl->leftLength * (maxCategories + 1) * sizeof(double));
      localTree->partitionData[model].right             = (double *)malloc_aligned((size_t)pl->rightLength * (maxCategories + 1) * sizeof(double));
      localTree->partitionData[model].EIGN              = (double*)malloc((size_t)pl->eignLength * sizeof(double));
      localTree->partitionData[model].EV                = (double*)malloc_aligned((size_t)pl->evLength * sizeof(double));
      localTree->partitionData[model].EI                = (double*)malloc((size_t)pl->eiLength * sizeof(double));

      localTree->partitionData[model].substRates        = (double *)malloc((size_t)pl->substRatesLength * sizeof(double));
      localTree->partitionData[model].frequencies       = (double*)malloc((size_t)pl->frequenciesLength * sizeof(double));
      localTree->partitionData[model].empiricalFrequencies       = (double*)malloc((size_t)pl->frequenciesLength * sizeof(double));
      localTree->partitionData[model].tipVector         = (double *)malloc_aligned((size_t)pl->tipVectorLength * sizeof(double));
      localTree->partitionData[model].symmetryVector    = (int *)malloc((size_t)pl->symmetryVectorLength  * sizeof(int));
      localTree->partitionData[model].frequencyGrouping = (int *)malloc((size_t)pl->frequencyGroupingLength  * sizeof(int));

      localTree->partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * maxCategories);

      localTree->partitionData[model].nonGTR = FALSE;            

      localTree->partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);      
      localTree->partitionData[model].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * ((size_t)localTree->mxtips + 1));


      localTree->partitionData[model].xVector = (double **)malloc(sizeof(double*) * (size_t)localTree->mxtips);   

      for(j = 0; j < (size_t)localTree->mxtips; j++)	        	  	  	  	 
	localTree->partitionData[model].xVector[j]   = (double*)NULL;   

      localTree->partitionData[model].xSpaceVector = (size_t *)calloc((size_t)localTree->mxtips, sizeof(size_t));  

      localTree->partitionData[model].sumBuffer = (double *)malloc_aligned(width *
									   (size_t)(localTree->partitionData[model].states) *
									   discreteRateCategories(localTree->rateHetModel) *
									   sizeof(double));


      /* data structure to store the marginal ancestral probabilities in the sequential version or for each thread */

      localTree->partitionData[model].ancestralBuffer = (double *)malloc_aligned(width *
										 (size_t)(localTree->partitionData[model].states) * 
										 sizeof(double));

      /* count and accumulate how many bytes we will need for storing a full ancestral vector. for this we addf over the per-partition space requirements in bytes */
      /* ancestralVectorWidth += ((size_t)(tr->partitionData[model].upper - tr->partitionData[model].lower) * (size_t)(localTree->partitionData[model].states) * sizeof(double)); */
      ancestralVectorWidth += ((size_t)(localTree->partitionData[model].upper - localTree->partitionData[model].lower) * (size_t)(localTree->partitionData[model].states) * sizeof(double));
      /* :TODO: do we have to use the original tree for that   */

      localTree->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      localTree->partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));

      if(width > 0 && localTree->saveMemory)
	{
	  localTree->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	  assert(4 == sizeof(unsigned int));
	  localTree->partitionData[model].gapVector = (unsigned int*)calloc((size_t)localTree->partitionData[model].gapVectorLength * 2 * (size_t)localTree->mxtips, sizeof(unsigned int));	  	    	  	  
	  localTree->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)localTree->mxtips) *								      
									       ((size_t)(localTree->partitionData[model].states)) *
									       discreteRateCategories(localTree->rateHetModel) * sizeof(double));
	}
      else
	{
	  localTree->partitionData[model].gapVectorLength = 0;  
	  localTree->partitionData[model].gapVector = (unsigned int*)NULL; 
	  localTree->partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }
}




