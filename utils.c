/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include "cycle.h"

#ifdef  _FINE_GRAIN_MPI
#include <mpi.h>
#endif



#ifdef _USE_PTHREADS
#include <pthread.h>

#endif

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
   special bug fix, enforces denormalized numbers to be flushed to zero,
   without this program is a tiny bit faster though.
#include <emmintrin.h> 
#define MM_DAZ_MASK    0x0040
#define MM_DAZ_ON    0x0040
#define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"
#include "globalVariables.h"


#define _PORTABLE_PTHREADS

boolean setupTree (tree *tr);

/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(tree *tr)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->td[0].executeModel[model] = tr->executeModel[model];
}

void storeValuesInTraversalDescriptor(tree *tr, double *value)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->td[0].parameterValues[model] = value[model];
}




static void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}


void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;

  int 
    res;


#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
     */

  ptr = malloc(size);

  if(ptr == (void*)NULL) 
    assert(0);

#ifdef __AVX
  assert(0);
#endif


#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 

  return ptr;
}





/* Marked for deletion --- Ask Alexi */

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
  {    
    case TREE_EVALUATION:
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, SUMMARIZE_LH, FALSE, FALSE);

      logFile = myfopen(temporaryFileName, "wb");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      if(adef->perGeneBranchLengths)
        printTreePerGene(tr, adef, temporaryFileName, "wb");
      break;
    case BIG_RAPID_MODE:     
      if(finalPrint)
      {
        switch(tr->rateHetModel)
        {
          case GAMMA:
          case GAMMA_I:
            Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint,
                SUMMARIZE_LH, FALSE, FALSE);

            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            if(adef->perGeneBranchLengths)
              printTreePerGene(tr, adef, temporaryFileName, "wb");
            break;
          case CAT:
            /*Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
              NO_BRANCHES, FALSE, FALSE);*/



            Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
                TRUE, SUMMARIZE_LH, FALSE, FALSE);




            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            break;
          default:
            assert(0);
        }
      }
      else
      {
        Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint,
            NO_BRANCHES, FALSE, FALSE);
        logFile = myfopen(temporaryFileName, "wb");
        fprintf(logFile, "%s", tr->tree_string);
        fclose(logFile);
      }    
      break;
    default:
      printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
  }
}

void printLog(tree *tr)
{
  FILE *logFile;
  double t;


  t = gettime() - masterTime;

  logFile = myfopen(logFileName, "ab");

  fprintf(logFile, "%f %f\n", t, tr->likelihood);

  fclose(logFile);


}




static void printBoth(FILE *f, const char* format, ... )
{
  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);
}

void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}




boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}



/* return the number of states of the given dataType,
   i.e. for DNA_DATA return 4 */
int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}


int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

const partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
              states    = p->states,
              tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  /*pLength.leftLength = pLength.rightLength = states * states;
    pLength.eignLength = states;
    pLength.evLength   = states * states;
    pLength.eiLength   = states * states;
    pLength.substRatesLength = (states * states - states) / 2;
    pLength.frequenciesLength = states;
    pLength.tipVectorLength   = tipLength * states;
    pLength.symmetryVectorLength = (states * states - states) / 2;
    pLength.frequencyGroupingLength = states;
    pLength.nonGTR = FALSE;*/

  return (&pLengths[dataType]); 
}

/* End of marked for deletion */







/* Number of evolutionary rate category used per site*/
size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
  {
    case CAT: /* PSR per site rates*/
      result = 1;
      break;
    case GAMMA:
      result = 4;
      break;
    default:
      assert(0);
  }

  return result;
}



double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static int filexists(char *filename)
{
  FILE 
    *fp = fopen(filename,"rb");

  int res; 

  if(fp)
  {
    res = 1;
    fclose(fp);
  }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
  {
    if(fp)
      return fp;
    else
    {	  
      printf("\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n", path);
      exit(-1);
      return (FILE *)NULL;
    }
  }
  else
  {
    if(fp)
      return fp;
    else
    {	 
      printf("\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
          path, mode);
      exit(-1);
      return (FILE *)NULL;
    }
  }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}









/* ensure p->x = 1 and p->next->x = 0 and p->next->next->x = 0 */
/* p->x marks the orientation of the node */
void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
  {
    p->x = s->x;
    s->x = 0;
  }

  assert(p->x);
}




/* connect p and q,  both share now the same branch length and p->back == q and p->back->back == p */
void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


/***********************reading and initializing input ******************/







boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
      y = 362436069,
      z = 21288629,
      w = 14921776,
      c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}
/* read the binary MSA format and store date on tree *tr */
void read_msa(tree *tr, const char *filename)
  {
    size_t 
      i,
      model;

    unsigned char *y;
  double **empiricalFrequencies;

    FILE 
      *byteFile = myfopen(filename, "rb");	 

    /* read the alignment info */
    myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
    myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
    myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
    myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);
    /* initialize topology */
    setupTree(tr);

    /* Joint branch length estimate is activated by default */
    /*
    if(adef->perGeneBranchLengths)
      tr->numBranches = tr->NumberOfModels;
    else
      tr->numBranches = 1;
    */
    tr->numBranches = 1;

    /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
       let's not worry about this right now, because it is indeed RAxML-specific */

    tr->aliaswgt                   = (int *)malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);	       

    tr->rateCategory    = (int *)    malloc((size_t)tr->originalCrunchedLength * sizeof(int));	  
    tr->wr              = (double *) malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
    tr->wr2             = (double *) malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
    tr->patrat          = (double*)  malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->patratStored    = (double*)  malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 
    tr->lhs             = (double*)  malloc((size_t)tr->originalCrunchedLength * sizeof(double)); 

    tr->executeModel   = (boolean *)malloc(sizeof(boolean) * (size_t)tr->NumberOfModels);

    for(i = 0; i < (size_t)tr->NumberOfModels; i++)
      tr->executeModel[i] = TRUE;



    /* data structures for convergence criterion need to be initialized after! setupTree */
    if(tr->searchConvergenceCriterion)
    {                     
      tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
      tr->h = initHashTable(tr->mxtips * 4);        
    }

    /* read tip names */
    for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int len;
      myBinFread(&len, sizeof(int), 1, byteFile);
      tr->nameList[i] = (char*)malloc(sizeof(char) * (size_t)len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      /*printf("%s \n", tr->nameList[i]);*/
    }  

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      addword(tr->nameList[i], tr->nameHash, i);

    /* read partition info (boudaries, data type) */
    empiricalFrequencies = (double **)malloc(sizeof(double *) * (size_t)tr->NumberOfModels);
    for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      int 
        len;

      pInfo 
        *p = &(tr->partitionData[model]);	   

      myBinFread(&(p->states),             sizeof(int), 1, byteFile);
      myBinFread(&(p->maxTipStates),       sizeof(int), 1, byteFile);
      myBinFread(&(p->lower),              sizeof(int), 1, byteFile);
      myBinFread(&(p->upper),              sizeof(int), 1, byteFile);
      myBinFread(&(p->width),              sizeof(int), 1, byteFile);
      myBinFread(&(p->dataType),           sizeof(int), 1, byteFile);
      myBinFread(&(p->protModels),         sizeof(int), 1, byteFile);
      myBinFread(&(p->autoProtModels),     sizeof(int), 1, byteFile);
      myBinFread(&(p->protFreqs),          sizeof(int), 1, byteFile);
      myBinFread(&(p->nonGTR),             sizeof(boolean), 1, byteFile);
      myBinFread(&(p->numberOfCategories), sizeof(int), 1, byteFile); 

      /* later on if adding secondary structure data

         int    *symmetryVector;
         int    *frequencyGrouping;
         */

      myBinFread(&len, sizeof(int), 1, byteFile);
      p->partitionName = (char*)malloc(sizeof(char) * (size_t)len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);

      empiricalFrequencies[model] = (double *)malloc(sizeof(double) * (size_t)tr->partitionData[model].states);
      myBinFread(empiricalFrequencies[model], sizeof(double), tr->partitionData[model].states, byteFile);	   
    }

    /* Read all characters from tips */
    y = (unsigned char *)malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

    tr->yVector = (unsigned char **)malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength];	

    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    fclose(byteFile);
  }
/* removed the static keyword for using this function in the examples */
boolean setupTree (tree *tr)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  size_t
    tips,
    inter; 


  tr->bigCutoff = FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tr->partitionContributions = (double *)malloc(sizeof(double) * (size_t)tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->partitionContributions[i] = -1.0;

  tr->perPartitionLH = (double *)malloc(sizeof(double) * (size_t)tr->NumberOfModels);


  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->perPartitionLH[i] = 0.0;	    



  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);




  tr->fracchanges  = (double *)malloc((size_t)tr->NumberOfModels * sizeof(double));




  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/



  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * (size_t)tr->NumberOfModels);
  tr->td[0].parameterValues = (double *)malloc(sizeof(double) * (size_t)tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;

  tr->constraintVector = (int *)malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;

    p->hash   =  KISS32(); /* hash table stuff */
    p->x      =  0;
    p->xBips  = 0;
    p->number =  i;
    p->next   =  p;
    p->back   = (node *)NULL;
    p->bInf   = (branchInfo *)NULL;            
    tr->nodep[i] = p;
  }

  for (i = tips + 1; i <= tips + inter; i++)
  {
    q = (node *) NULL;
    for (j = 1; j <= 3; j++)
    {	 
      p = p0++;
      if(j == 1)
      {
        p->xBips = 1;
        p->x = 1;
      }
      else
      {
        p->xBips = 0;
        p->x =  0;
      }
      p->number = i;
      p->next   = q;
      p->bInf   = (branchInfo *)NULL;
      p->back   = (node *) NULL;
      p->hash   = 0;       
      q = p;
    }
    p->next->next->next = p;
    tr->nodep[i] = p;
  }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;

  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  tr->partitionData = (pInfo*)malloc(sizeof(pInfo) * (size_t)tr->NumberOfModels);

  return TRUE;
}













/* Marked for deletion */

static void initAdef(analdef *adef)
{   
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = 10;
  adef->bestTrav               = 10;
  adef->initialSet             = FALSE; 
  adef->mode                   = BIG_RAPID_MODE; 
  adef->likelihoodEpsilon      = 0.1;

  adef->permuteTreeoptimize    = FALSE; 
  adef->perGeneBranchLengths   = FALSE;  

  adef->useCheckpoint          = FALSE;

#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
  adef->num_generations        = 10000;
#endif
}

/* end of marked for deletion */



static boolean modelExists(char *model, tree *tr)
{
  /********** BINARY ********************/

  if(strcmp(model, "PSR") == 0)
  {
    tr->rateHetModel = CAT;
    return TRUE;
  }

  if(strcmp(model, "GAMMA") == 0)
  {
    tr->rateHetModel = GAMMA;
    return TRUE;
  }


  return FALSE;
}


/* Marked for deletion */
static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
  {
    if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
      return -1;
  }
  else
  {
    if(strcmp(argv[*optind], "--") == 0)
    {
      *optind =  *optind + 1;
      return -1;
    }
  }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
  {
    printf("\n Error: illegal option \"-%c\" \n\n", c);
    if(argv[*optind][++sp] == '\0')
    {
      *optind =  *optind + 1;
      sp = 1;
    }
    return('?');
  }
  if(*++cp == ':')
  {
    if(argv[*optind][sp+1] != '\0')
    {
      *optarg = &argv[*optind][sp+1];
      *optind =  *optind + 1;
    }
    else
    {
      *optind =  *optind + 1;
      if(*optind >= argc)
      {
        printf("\n Error: option \"-%c\" requires an argument\n\n", c);
        sp = 1;
        return('?');
      }
      else
      {
        *optarg = argv[*optind];
        *optind =  *optind + 1;
      }
    }
    sp = 1;
  }
  else
  {
    if(argv[*optind][++sp] == '\0')
    {
      sp = 1;
      *optind =  *optind + 1;
    }
    *optarg = 0;
  }
  return(c);
}
/* End of marked for deletion */


void init_defaults (tree * tr)
{
  /*********** tr inits **************/

  tr->numberOfThreads = 1; 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;

  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;

  tr->startingTree = randomTree;
  tr->randomNumberSeed = 12345;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 
  tr->useMedian = FALSE;
  /* recom */
  tr->useRecom = FALSE;
  tr->rvec = (recompVectors*)NULL;
  /* recom */

  /********* tr inits end*************/

}


/*********************************** *********************************************************/












/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/



















void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
  {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
  }
}







/************************************************************************************/


#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))




boolean isThisMyPartition(tree *localTree, int tid, int model)
{ 
  if(localTree->partitionAssignment[model] == tid)
    return TRUE;
  else
    return FALSE;
}

static void computeFractionMany(tree *localTree, int tid)
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



static void computeFraction(tree *localTree, int tid, int n)
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












#ifdef _USE_PTHREADS



inline static void broadcastTraversalInfo(tree *localTree, tree *tr)
{
  /* the one below is a hack we are re-assigning the local pointer to the global one
     the memcpy version below is just for testing and preparing the
     fine-grained MPI BlueGene version */

  /* TODO: we should reset this at some point, the excplicit copy is just done for testing */

  if(0)
    {
      localTree->td[0] = tr->td[0];
    }
  else
    {
      localTree->td[0].count               = tr->td[0].count;
      localTree->td[0].functionType        = tr->td[0].functionType;
      localTree->td[0].traversalHasChanged = tr->td[0].traversalHasChanged;
     
      /* memcpy -> memmove (see ticket #43). This function is sometimes called with localTree == tr,
       * in which case some memcpy implementations can corrupt the buffers.
       */
      
      memmove(localTree->td[0].executeModel,    tr->td[0].executeModel,    sizeof(boolean) * localTree->NumberOfModels);
      memmove(localTree->td[0].parameterValues, tr->td[0].parameterValues, sizeof(double) * localTree->NumberOfModels);
      
      if(localTree->td[0].traversalHasChanged)
	memmove(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));
    }
}


static void collectDouble(double *dst, double *src, tree *tr, int n, int tid)
{
  int 
    model,
    i;

  if(tr->manyPartitions)
    for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(isThisMyPartition(tr, tid, model))	
        for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
          dst[i] = src[i];       
    }
  else
    for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
      {
        if(i % n == tid)
          dst[i] = src[i];
      }
    }
}


static void broadcastPerSiteRates(tree *tr, tree *localTree)
{
  int
    i = 0,
      model = 0;  

  for(model = 0; model < localTree->NumberOfModels; model++)
  {
    localTree->partitionData[model].numberOfCategories = tr->partitionData[model].numberOfCategories;

    for(i = 0; i < localTree->partitionData[model].numberOfCategories; i++)
      localTree->partitionData[model].perSiteRates[i] = tr->partitionData[model].perSiteRates[i];
  }

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


  for(model = 0; model < localTree->NumberOfModels; model++)
    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
}

static void broadCastRates(tree *localTree, tree *tr, int tid)
{
  int 
    model;
#ifdef _LOCAL_DISCRETIZATION        
  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
      
      if(tid > 0)
	memcpy(localTree->partitionData[model].substRates,        tr->partitionData[model].substRates, pl->substRatesLength * sizeof(double));
      
      initReversibleGTR(localTree, model);     
    }    
#else
  if(tid > 0)
  {	 
    for(model = 0; model < localTree->NumberOfModels; model++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));

      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));		  
      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));	      	      	     
    }
  }
#endif
}



static void broadCastAlpha(tree *localTree, tree *tr, int tid)
{
  int 
    model; 

#ifdef _LOCAL_DISCRETIZATION 
  for(model = 0; model < localTree->NumberOfModels; model++)
  {
    localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
    makeGammaCats(localTree->partitionData[model].alpha, localTree->partitionData[model].gammaRates, 4, tr->useMedian); 
  }   
#else
  if(tid > 0)
  {	 
    for(model = 0; model < localTree->NumberOfModels; model++)
      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
  }
#endif	      	      	     
}

static void initializePartitions(tree *tr, tree *localTree, int tid, int n);
/* this function here handles all parallel regions in the Pthreads version, when we enter 
   this function masterBarrier() has ben called by the master thread from within the sequential 
   part of the program, tr is the tree at the master thread, localTree the tree at the worker threads

   While this is not necessary, adress spaces of threads are indeed separated for easier transition to 
   a distributed memory paradigm 
   */



static void execFunction(tree *tr, tree *localTree, int tid, int n)
{
  int
    i,
    currentJob,   
    model,
    localCounter;

  /* some stuff associated with the barrier implementation using Pthreads and busy wait */

  currentJob = threadJob >> 16;

  /* here the master sends and all threads/processes receive the traversal descriptor */

  broadcastTraversalInfo(localTree, tr);

  /* make sure that nothing is going wrong */

  assert(currentJob == localTree->td[0].functionType);

  /*  switch over the function type in the traversal descriptor to figure out which 
      parallel region to execute */

  switch(localTree->td[0].functionType)
  {            
    case THREAD_NEWVIEW:      
      /* just a newview on the fraction of sites that have been assigned to this thread */

      newviewIterative(localTree, 0);
      break;     
    case THREAD_EVALUATE:            

      reduceEvaluateIterative(localTree, tid);            	        
      break;	
    case THREAD_MAKENEWZ_FIRST:
    case  THREAD_MAKENEWZ:
      /* this is the first call from within makenewz that requires getting the likelihood vectors to the left and 
         right of the branch via newview and doing som eprecomputations.

         For details see comments in makenewzGenericSpecial.c 
         */

      {
        int 
          b;

        volatile double
          dlnLdlz[NUM_BRANCHES],
          d2lnLdlz2[NUM_BRANCHES];	

        if(localTree->td[0].functionType == THREAD_MAKENEWZ_FIRST)
          makenewzIterative(localTree);	
        execCore(localTree, dlnLdlz, d2lnLdlz2);

        /* gather the first and second derivatives that have been written by each thread */
        /* as for evaluate above, the final sum over the derivatives will be computed by the 
           master thread in its sequential part of the code */


        for(b = 0; b < localTree->numBranches; b++)
        {	     
          reductionBuffer[tid * localTree->numBranches + b]    = dlnLdlz[b];
          reductionBufferTwo[tid * localTree->numBranches + b] = d2lnLdlz2[b];
        }	
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
      broadCastAlpha(localTree, tr, tid);

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

      /* need to copy base freqs before local per-thread/per-process Q matrix exponentiation */
#ifdef _LOCAL_DISCRETIZATION   
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
      broadCastAlpha(localTree, tr, tid); 

      /*
         copy initial model parameters, the Q matrix and alpha are initially, when we start our likelihood search 
         set to default values. 
         Hence we need to copy all those values that are required for computing the likelihood 
         with newview(), evaluate() and makenez() to the private memory of the threads 
         */


      if(tid > 0 && localTree->rateHetModel == CAT)
      {	  	  
        for(model = 0; model < localTree->NumberOfModels; model++)	    	     
          localTree->partitionData[model].numberOfCategories      = tr->partitionData[model].numberOfCategories;	    

        /* this is only relevant for the PSR model, we can worry about this later */

        memcpy(localTree->patrat,       tr->patrat,      localTree->originalCrunchedLength * sizeof(double));
        memcpy(localTree->patratStored, tr->patratStored, localTree->originalCrunchedLength * sizeof(double));	  
      }          

      break;    
    case THREAD_RATE_CATS:            
      /* this is for optimizing per-site rate categories under PSR, let's worry about this later */

      if(tid > 0)
      {
        localTree->lower_spacing = tr->lower_spacing;
        localTree->upper_spacing = tr->upper_spacing;
      }

      optRateCatPthreads(localTree, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

      if(tid > 0)
      {
        collectDouble(tr->patrat,       localTree->patrat,         localTree, n, tid);
        collectDouble(tr->patratStored, localTree->patratStored,   localTree, n, tid);
        collectDouble(tr->lhs,                localTree->lhs,                  localTree, n, tid);
      }
      break;
    case THREAD_COPY_RATE_CATS:
      /* this is invoked when we have changed the per-site rate category assignment
         In essence it distributes the new per site rates to all threads 
         */

      if(tid > 0)
      {	  
        memcpy(localTree->patrat,       tr->patrat,         localTree->originalCrunchedLength * sizeof(double));
        memcpy(localTree->patratStored, tr->patratStored,   localTree->originalCrunchedLength * sizeof(double));
        broadcastPerSiteRates(tr, localTree);
      }

      for(model = 0; model < localTree->NumberOfModels; model++)
      {
        localTree->partitionData[model].numberOfCategories = tr->partitionData[model].numberOfCategories;

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
      break;
    case THREAD_PER_SITE_LIKELIHOODS:

      /* compute per-site log likelihoods for the sites/partitions 
         that are handled by this thread */

      perSiteLogLikelihoodsPthreads(localTree, localTree->lhs, n, tid);

      /* do a parallel gather operation, the threads will write their results 
         into the global buffer tr->lhs that will then contain all per-site log likelihoods
         in the proper order 
         */

      if(tid > 0)
        collectDouble(tr->lhs,                localTree->lhs,                  localTree, n, tid);

      break;
      /* check for errors */

  case THREAD_NEWVIEW_ANCESTRAL:
    /* This one here is easy, since we just need to invoke the function for calculating ancestral states here */
    newviewAncestralIterative(localTree);
    break;
  case THREAD_GATHER_ANCESTRAL:
      {
	/* gather the ancestral states from the threads */

	double
	  *contigousVector = tr->ancestralVector;
	
	double 
	  *stridedVector = localTree->partitionData[model].ancestralBuffer;
	
	size_t	 
	  globalColumnCount = 0,
	  globalCount       = 0;	
	
	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    size_t
	      blockRequirements = (size_t)(localTree->partitionData[model].states);
	    
	    if(localTree->manyPartitions)
	      {	   
		size_t 
		  width = (localTree->partitionData[model].upper - localTree->partitionData[model].lower) * blockRequirements;
		
		if(isThisMyPartition(localTree, tid, model))		 
		  memcpy(&contigousVector[globalCount], stridedVector, sizeof(double) * width);			 

		globalCount += width;
	      }
	    else
	      {				
		size_t	     	     	   
		  localCount = 0;	   
					    	   	    		
		for(globalColumnCount = localTree->partitionData[model].lower; globalColumnCount < localTree->partitionData[model].upper; globalColumnCount++)
		  {	
		    
		    if(globalColumnCount % n == tid)
		      {			
			memcpy(&contigousVector[globalCount], &stridedVector[localCount], sizeof(double) * blockRequirements);		
			
			localCount += blockRequirements;
		      }	
		    
		    globalCount += blockRequirements;
		  }	    
		
		assert(localCount == (localTree->partitionData[model].width * (int)blockRequirements));
	      }
	  }	
      }
    break;
  default:
      printf("Job %d\n", currentJob);
      assert(0);
  }
}




void masterBarrier(int jobType, tree *tr)
{
  const int 
    n = tr->numberOfThreads;

  int 
    i, 
    sum;

  tr->td[0].functionType = jobType;

  jobCycle = !jobCycle;
  threadJob = (jobType << 16) + jobCycle;

  execFunction(tr, tr, 0, n);


  do
  {
    for(i = 1, sum = 1; i < n; i++)
      sum += barrierBuffer[i];
  }
  while(sum < n);  

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
}

#ifndef _PORTABLE_PTHREADS

static void pinToCore(int tid)
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

static void *likelihoodThread(void *tData)
{
  threadData *td = (threadData*)tData;
  tree
    *tr = td->tr,
    *localTree = (tree *)malloc(sizeof(tree));
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

  return (void*)NULL;
}

static void startPthreads(tree *tr)
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

  reductionBuffer          = (volatile double *)malloc(sizeof(volatile double) *  (size_t)tr->numberOfThreads * (size_t)tr->NumberOfModels);
  reductionBufferTwo       = (volatile double *)malloc(sizeof(volatile double) *  (size_t)tr->numberOfThreads * (size_t)tr->NumberOfModels);  
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


#endif











#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))    

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

static void multiprocessorScheduling(tree *tr, int tid)
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

#endif

static void initializePartitions(tree *tr, tree *localTree, int tid, int n)
{ 
  size_t
    /* in ancestralVectorWidth we store the total length in bytes (!) of 
       one conditional likelihood array !
       we need to know this length such that in the pthreads version the master thread can actually 
       gather the scattered ancestral probabilities from the threads such that they can be printed to screen!
    */
       
    ancestralVectorWidth = 0,
    model,
    maxCategories;

  localTree->threadID = tid; 


#ifdef _USE_PTHREADS
  if(tid > 0)
  {
    int totalLength = 0;

    assert(localTree != tr);

    localTree->numberOfThreads         = tr->numberOfThreads;
    localTree->manyPartitions          = tr->manyPartitions;
    localTree->NumberOfModels          = tr->NumberOfModels;           
    localTree->rateHetModel            = tr->rateHetModel;
    localTree->useMedian               = tr->useMedian;
    localTree->saveMemory              = tr->saveMemory;
    localTree->maxCategories           = tr->maxCategories;      
    localTree->originalCrunchedLength  = tr->originalCrunchedLength;    
    localTree->mxtips                  = tr->mxtips;     
    localTree->numBranches             = tr->numBranches;
#ifdef _LOCAL_DISCRETIZATION
    localTree->aliaswgt                = tr->aliaswgt;
#endif    
    localTree->lhs                     = (double*)malloc(sizeof(double)   * (size_t)localTree->originalCrunchedLength);     
    localTree->perPartitionLH          = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);     
    localTree->fracchanges             = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);
    localTree->partitionContributions  = (double*)malloc(sizeof(double)   * (size_t)localTree->NumberOfModels);
    localTree->partitionData           = (pInfo*)malloc(sizeof(pInfo)     * (size_t)localTree->NumberOfModels);

    localTree->td[0].count = 0;
    localTree->td[0].ti              = (traversalInfo *)malloc(sizeof(traversalInfo) * (size_t)localTree->mxtips);
    localTree->td[0].executeModel    = (boolean *)malloc(sizeof(boolean) * (size_t)localTree->NumberOfModels);
    localTree->td[0].parameterValues = (double *)malloc(sizeof(double) * (size_t)localTree->NumberOfModels);

    localTree->patrat       = (double*)malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);
    localTree->patratStored = (double*)malloc(sizeof(double) * (size_t)localTree->originalCrunchedLength);      

    for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      localTree->partitionData[model].numberOfCategories    = tr->partitionData[model].numberOfCategories;
      localTree->partitionData[model].states                = tr->partitionData[model].states;
      localTree->partitionData[model].maxTipStates          = tr->partitionData[model].maxTipStates;
      localTree->partitionData[model].dataType              = tr->partitionData[model].dataType;
      localTree->partitionData[model].protModels            = tr->partitionData[model].protModels;
      localTree->partitionData[model].protFreqs             = tr->partitionData[model].protFreqs;	  
      localTree->partitionData[model].lower                 = tr->partitionData[model].lower;
      localTree->partitionData[model].upper                 = tr->partitionData[model].upper;	 
      localTree->perPartitionLH[model]                      = 0.0;

      totalLength += (localTree->partitionData[model].upper -  localTree->partitionData[model].lower);
    }

    assert(totalLength == localTree->originalCrunchedLength);
    /* recomp */
    localTree->useRecom = tr->useRecom;
    localTree->vectorRecomFraction = tr->vectorRecomFraction;
    /* E recomp */
  }

  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    localTree->partitionData[model].width        = 0;

  if(tr->manyPartitions)    
    multiprocessorScheduling(localTree, tid);

  if(tr->manyPartitions)
    computeFractionMany(localTree, tid);
  else
    computeFraction(localTree, tid, n);



#else
  assert(tr == localTree);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    assert(tr->partitionData[model].width == tr->partitionData[model].upper - tr->partitionData[model].lower);
#endif	   

  maxCategories = (size_t)localTree->maxCategories;

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

    ancestralVectorWidth += ((size_t)(tr->partitionData[model].upper - tr->partitionData[model].lower) * (size_t)(localTree->partitionData[model].states) * sizeof(double));
			     
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

#ifdef _USE_PTHREADS
  /* we only need to allocate this buffer for gathering the data at the master thread */
      if(tid == 0)
	tr->ancestralVector = (double *)malloc(ancestralVectorWidth);
#endif



#if ! (defined(_USE_PTHREADS) || defined(_FINE_GRAIN_MPI))
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
      offset,
      countOffset,
      myLength = 0;

    for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
      myLength += localTree->partitionData[model].width;         

    /* assign local memory for storing sequence data */

    localTree->y_ptr = (unsigned char *)malloc(myLength * (size_t)(localTree->mxtips) * sizeof(unsigned char));
    assert(localTree->y_ptr != NULL);

    for(i = 0; i < (size_t)localTree->mxtips; i++)
    {
      for(model = 0, offset = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
      {
        localTree->partitionData[model].yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
        countOffset +=  localTree->partitionData[model].width;
      }
      assert(countOffset == myLength);
    }

    /* figure in data */

    if(tr->manyPartitions)
      for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
      {
        if(isThisMyPartition(localTree, tid, model))
        {
          assert(localTree->partitionData[model].upper - localTree->partitionData[model].lower == localTree->partitionData[model].width);

          for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++, localCounter++)
          {	    
            localTree->partitionData[model].wgt[localCounter]          = tr->aliaswgt[globalCounter];	      	     

            for(j = 1; j <= (size_t)localTree->mxtips; j++)
              localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter]; 	     

            globalCounter++;
          }
        }
        else
          globalCounter += (localTree->partitionData[model].upper - localTree->partitionData[model].lower);
      }
    else
      for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
      {
        for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++)
        {
          if(i % (size_t)n == (size_t)tid)
          {
            localTree->partitionData[model].wgt[localCounter]          = tr->aliaswgt[globalCounter];	      	     		 

            for(j = 1; j <= (size_t)localTree->mxtips; j++)
              localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter]; 	     

            localCounter++;
          }
          globalCounter++;
        }
      }    
  }
#endif

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
}

/* Pick up a random inner subtree */
nodeptr pickRandomSubtree(tree *tr)
{
  nodeptr p;
  do
  {
    int exitDirection = rand() % 3; 
    p = tr->nodep[(rand() % (tr->mxtips - 2)) + 1 + tr->mxtips];
    switch(exitDirection)
    {
      case 0:
        break;
      case 1:
        p = p->next;
        break;
      case 2:
        p = p->next->next;
        break;
      default:
        assert(0);
    }
  }
  while(isTip(p->next->back->number, tr->mxtips) && isTip(p->next->next->back->number, tr->mxtips));
  assert(!isTip(p->number, tr->mxtips));
  return p;
}

/* small example program that executes ancestral state computations 
   on the entire subtree rooted at p.

   Note that this is a post-order traversal.
*/

  
static void computeAllAncestralVectors(nodeptr p, tree *tr)
{
  /* if this is not a tip, for which evidently it does not make sense 
     to compute the ancestral sequence because we have the real one ....
  */

  if(!isTip(p->number, tr->mxtips))
    {
      /* descend recursively to compute the ancestral states in the left and right subtrees */

      computeAllAncestralVectors(p->next->back, tr);
      computeAllAncestralVectors(p->next->next->back, tr);
      
      /* then compute the ancestral state at node p */

      newviewGenericAncestral(tr, p);

      /* and print it to terminal, the two booleans that are set to true here 
	 tell the function to print the marginal probabilities as well as 
	 a discrete inner sequence, that is, ACGT etc., always selecting and printing 
	 the state that has the highest probability */

      printAncestralState(p, TRUE, TRUE, tr);
    }
}
