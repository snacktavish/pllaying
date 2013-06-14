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

/** @file utils.c
 *  
 *  @brief Miscellaneous general utility and helper functions
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
#include <assert.h>
#include "cycle.h"
#include "parser/phylip/ssort.h"


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
//#include "phylip_parser/lexer.h"
//#include "phylip_parser/phylip.h"
//#include "phylip_parser/xalloc.h"
//#include "phylip_parser/msa_sites.h"


#include "globalVariables.h"
#include "mem_alloc.h"
#include "queue.h"
#include "parser/partition/part.h"
#include "parser/phylip/phylip.h"
#include "parser/newick/newick.h"
#include "utils.h"


extern unsigned int mask32[32];


/***************** UTILITY FUNCTIONS **************************/

void storeExecuteMaskInTraversalDescriptor(pllInstance *tr, partitionList *pr)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].executeModel[model] = pr->partitionData[model]->executeModel;

}

void storeValuesInTraversalDescriptor(pllInstance *tr, partitionList *pr, double *value)
{
  int model;

  for(model = 0; model < pr->numberOfPartitions; model++)
    tr->td[0].parameterValues[model] = value[model];
}

#ifdef EXPERIMENTAL
void read_phylip_msa(pllInstance * tr, const char * filename, int format, int type)
{
    size_t
      i, j,
      model;

  struct phylip_data * pd;
  struct msa_sites * ms;
  double **empiricalFrequencies;

  pd = pl_phylip_parse (filename, format);

  ms = construct_msa_sites (pd, SITES_CREATE | SITES_COMPUTE_WEIGHTS);

  free_phylip_struct (pd);
  pd = transpose (ms);
  free_sites_struct (ms);



  tr->mxtips                 = pd->taxa;
  tr->originalCrunchedLength = pd->seqlen;
  pr->numberOfPartitions         = 1;

  setupTree(tr, PLL_TRUE);

  tr->gapyness               = 0.03;   /* number of undetermined chars / alignment size */

  /* TODO: The next two lines were commented in model-sep branch */
  tr->aliaswgt = pl_phylip_deldups (&pd);
  tr->originalCrunchedLength = pd->seqlen;
  pr->perGeneBranchLengths = PLL_FALSE;

  pl_phylip_subst (pd, DNA_DATA);          /* TODO: Change to reflect the input type */

  tr->rateCategory           =  (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));

  tr->patrat                 =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));
  tr->patratStored           =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));
  tr->lhs                    =  (double *) rax_malloc ((size_t)tr->originalCrunchedLength * sizeof (double));

  tr->executeModel   = (boolean *)rax_malloc(sizeof(boolean) * (size_t)pr->numberOfPartitions);



        
  for(i = 0; i < (size_t)pr->numberOfPartitions; i++)
    tr->executeModel[i] = PLL_TRUE;



  /* data structures for convergence criterion need to be initialized after! setupTree */
  if(tr->searchConvergenceCriterion)
  {
    tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
    tr->h = initHashTable(tr->mxtips * 4);
  }

  /* read tip names */
  for(i = 1; i <= (size_t)tr->mxtips; i++)
  {
    tr->nameList[i] = pd->label[i - 1];
  }

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  /* read partition info (boudaries, data type) */
  empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * (size_t)pr->numberOfPartitions);
  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
  {
    int
      len;

    pInfo
      *p = &(pr->partitionData[model]);

    p->states             =  4;   /* TODO: according to the type */
    p->maxTipStates       = 16;   /* TODO: according to the type */
    p->lower              =  0;
    p->upper              = pd->seqlen;
    p->width              = p->upper - p->lower;
    p->dataType           =   DNA_DATA; /* TODO: dna type */
    p->protModels         =   2;
    p->autoProtModels     =   0;
    p->protFreqs          =   0;
    p->nonGTR             =   PLL_FALSE;
    p->numberOfCategories =   0;
    
    /* later on if adding secondary structure data

       int    *symmetryVector;
       int    *frequencyGrouping;
       */

    p->partitionName = strdup ("PartName");

//    empiricalFrequencies[model] = (double *)malloc(sizeof(double) * (size_t)pr->partitionData[model]->states);
//    empiricalfrequencies[model][0] = 0.2036082474;
//    empiricalfrequencies[model][1] = 0.2268041237;
//    empiricalfrequencies[model][2] = 0.2731958763;
//    empiricalfrequencies[model][3] = 0.2963917526;
  }
  /* Read all characters from tips */
//  y = (unsigned char *)malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

  tr->yVector = (char **) rax_malloc(sizeof(char*) * (tr->mxtips+1));
 for (i=0; i < tr->mxtips; ++i)
        tr->yVector[i+1] = pd->seq[i]; //(unsigned char **)malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));
 
 #ifndef _USE_PTHREADS
 #ifndef _FINE_GRAIN_MPI
  //initializePartitionsSequential(tr); 
  initializePartitions (tr, tr, 0, 0);
 #endif
 #endif
}
#endif

/** @brief Read MSA from a file and setup the tree
 *
 *  Reads the MSA from \a filename and constructs
 *  the tree \a tr and sets up partition and model data
 *
 *  @todo This will be soon replaced by \a read_phylip_msa
 *
 *  @param tr
 *    Pointer to the tree instance to be set up
 *
 *  @param filename
 *    Filename containing the MSA
 *
 */
void read_msa(pllInstance *tr, partitionList *pr, const char *filename)
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
    myBinFread(&(pr->numberOfPartitions),         sizeof(int), 1, byteFile);

    /* initialize topology */

    /* Joint branch length estimate is activated by default */
    /*
    if(adef->perGeneBranchLengths)
      tr->numBranches = pr->numberOfPartitions;
    else
      tr->numBranches = 1;
    */
    pr->perGeneBranchLengths = PLL_FALSE;
    setupTree(tr, PLL_TRUE, pr);
    
    myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);

    /* If we use the RF-based convergence criterion we will need to allocate some hash tables.
       let's not worry about this right now, because it is indeed RAxML-specific */

    tr->aliaswgt                   = (int *)rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    myBinFread(tr->aliaswgt, sizeof(int), tr->originalCrunchedLength, byteFile);

    tr->rateCategory    = (int *)    rax_malloc((size_t)tr->originalCrunchedLength * sizeof(int));
    tr->patrat          = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->patratStored    = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
    tr->lhs             = (double*)  rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));

    for(i = 0; i < (size_t)pr->numberOfPartitions; i++)
      pr->partitionData[i]->executeModel = PLL_TRUE;



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
      tr->nameList[i] = (char*)rax_malloc(sizeof(char) * (size_t)len);
      myBinFread(tr->nameList[i], sizeof(char), len, byteFile);
      /*printf("%s \n", tr->nameList[i]);*/
    }

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      addword(tr->nameList[i], tr->nameHash, i);

    /* read partition info (boudaries, data type) */
    empiricalFrequencies = (double **)rax_malloc(sizeof(double *) * (size_t)pr->numberOfPartitions);
    for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    {
      int
        len;

      pInfo
        *p = pr->partitionData[model];

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
      p->partitionName = (char*)rax_malloc(sizeof(char) * (size_t)len);
      myBinFread(p->partitionName, sizeof(char), len, byteFile);

      empiricalFrequencies[model] = (double *)rax_malloc(sizeof(double) * (size_t)pr->partitionData[model]->states);
      myBinFread(empiricalFrequencies[model], sizeof(double), pr->partitionData[model]->states, byteFile);
    }
    /* Read all characters from tips */
    y = (unsigned char *)rax_malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));

    tr->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));

    for(i = 1; i <= (size_t)tr->mxtips; i++)
      tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength];

    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    /* Initialize the model */
    //printf("Here 1!\n");
    initializePartitionsSequential(tr, pr);
    //printf("Here 2!\n");
    initModel(tr, empiricalFrequencies, pr);
    fclose(byteFile);
  }



void myBinFread(void *ptr, size_t size, size_t nmemb, FILE *byteFile)
{  
  size_t
    bytes_read;

  bytes_read = fread(ptr, size, nmemb, byteFile);

  assert(bytes_read == nmemb);
}

#if 0
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

  ptr = rax_malloc(size);

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
#endif



/* Marked for deletion 
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

*/


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

void printResult(pllInstance *tr, partitionList *pr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
  {    
    case TREE_EVALUATION:
      Tree2String(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

      logFile = myfopen(temporaryFileName, "wb");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      if(adef->perGeneBranchLengths)
        printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
      break;
    case BIG_RAPID_MODE:
      if(finalPrint)
      {
        switch(tr->rateHetModel)
        {
          case GAMMA:
          case GAMMA_I:

            Tree2String(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint,
                PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);

            logFile = myfopen(temporaryFileName, "wb");
            fprintf(logFile, "%s", tr->tree_string);
            fclose(logFile);

            if(adef->perGeneBranchLengths)
              printTreePerGene(tr, pr, adef, temporaryFileName, "wb");
            break;
          case CAT:
            /*Tree2String(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint, adef,
              PLL_NO_BRANCHES, PLL_FALSE, PLL_FALSE);*/


            Tree2String(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE, PLL_TRUE, PLL_FALSE, PLL_FALSE,
                PLL_TRUE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);




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
        Tree2String(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE, PLL_TRUE, PLL_FALSE, PLL_FALSE, finalPrint,
            PLL_NO_BRANCHES, PLL_FALSE, PLL_FALSE);
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



/* Marked for deletion 
boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}
*/

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}

/*
int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}
*/

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}

/*
char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}
*/

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
    pLength.nonGTR = PLL_FALSE;*/
  return (&pLengths[dataType]); 
}

size_t discreteRateCategories(int rateHetModel)
{
  size_t 
    result;

  switch(rateHetModel)
  {
    case CAT:
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

/* Marked for deletion 
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
*/


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


/** @brief Check whether a node is a tip.
  * 
  * @param number
  *  Node number to be checked
  *
  * @param maxTips
  *  Number of tips in the tree
  *
  * @return
  *   \b PLL_TRUE if tip, \b PLL_FALSE otherwise
  */
boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return PLL_TRUE;
  else
    return PLL_FALSE;
}

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


/** @brief Connect two nodes and assign branch lengths 
  * 
  * Connect the two nodes \a p and \a q in each partition \e i with a branch of
  * length \a z[i]
  *
  * @param p
  *   Node \a p
  * 
  * @param q
  *   Node \a q
  *
  * @param numBranches
  *   Number of partitions
  */
void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

/* connects node p with q and assigns the branch lengths z for the whole vector*/
void hookupFull (nodeptr p, nodeptr q, double *z)
{
  //int i;

  p->back = q;
  q->back = p;

  memcpy(p->z, z, NUM_BRANCHES*sizeof(double) );
  memcpy(q->z, z, NUM_BRANCHES*sizeof(double) );
  //for(i = 0; i < numBranches; i++)
  //  p->z[i] = q->z[i] = z[i];

}

/* connect node p with q and assign the default branch lengths */
void hookupDefault (nodeptr p, nodeptr q)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < NUM_BRANCHES; i++)
    p->z[i] = q->z[i] = PLL_DEFAULTZ;

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

/* removed the static keyword for using this function in the examples */
boolean setupTree (pllInstance *tr, boolean doInit, partitionList *partitions)
{
  nodeptr  p0, p, q;
  int
    i,
    j;

  int
    tips,
    inter; 

  if(doInit)
    init_default(tr);

  tr->bigCutoff = PLL_FALSE;

  tr->maxCategories = MAX(4, tr->categories);

  tips  = (size_t)tr->mxtips;
  inter = (size_t)(tr->mxtips - 1);

  tr->treeStringLength = tr->mxtips * (PLL_NMLNGTH + 128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));


  /*TODO, must that be so long ?*/

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].executeModel = (boolean *)rax_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].parameterValues = (double *)rax_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  tr->fracchange = -1.0;

  tr->constraintVector = (int *)rax_malloc((2 * (size_t)tr->mxtips) * sizeof(int));

  tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));


  p0 = (nodeptr)rax_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr->nodeBaseAddress = p0;


  tr->nodep = (nodeptr *) rax_malloc((2* (size_t)tr->mxtips) * sizeof(nodeptr));
  assert(tr->nodep);    

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
  {
    p = p0++;

    p->hash   =  KISS32(); /* hast table stuff */
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

  tr->likelihood  = PLL_UNLIKELY;
  tr->start       = (node *) NULL;  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionSmoothed[i] = PLL_FALSE;

  tr->bitVectors = (unsigned int **)NULL;

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;

  tr->nameHash = initStringHashTable(10 * tr->mxtips);

  for (i = 0; i < partitions->numberOfPartitions; i++) {
	partitions->partitionData[i] = (pInfo*)rax_malloc (sizeof(pInfo));
	partitions->partitionData[i]->partitionContribution = -1.0;
	partitions->partitionData[i]->partitionLH = 0.0;
	partitions->partitionData[i]->fracchange = 1.0;
  }

  return PLL_TRUE;
}


boolean modelExists(char *model, pllInstance *tr)
{
  /********** BINARY ********************/

  if(strcmp(model, "PSR") == 0)
  {
    tr->rateHetModel = CAT;
    return PLL_TRUE;
  }

  if(strcmp(model, "GAMMA") == 0)
  {
    tr->rateHetModel = GAMMA;
    return PLL_TRUE;
  }


  return PLL_FALSE;
}



/*********************************** *********************************************************/


void init_default(pllInstance *tr)
{

  /*********** tr inits **************/

  tr->numberOfThreads = 1; 
  tr->doCutoff = PLL_TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = PLL_FALSE;
  tr->rateHetModel = GAMMA;

  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = PLL_FALSE;

  tr->fastScaling = PLL_FALSE;

  tr->manyPartitions = PLL_FALSE;

  tr->startingTree = randomTree;

  tr->categories             = 25;

  tr->grouped = PLL_FALSE;
  tr->constrained = PLL_FALSE;

  tr->gapyness               = 0.0; 
  tr->useMedian = PLL_FALSE;
  /* recom */
  tr->useRecom = PLL_FALSE;
  tr->rvec = (recompVectors*)NULL;
  /* recom */

  /********* tr inits end*************/

}








/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/



/* Delete it at some point */
void printLog(pllInstance *tr)
{
  FILE *logFile;
  double t;


  t = gettime() - masterTime;

  logFile = myfopen(logFileName, "ab");

  fprintf(logFile, "%f %f\n", t, tr->likelihood);

  fclose(logFile);


}


void getDataTypeString(pllInstance *tr, pInfo *partitionInfo, char typeOfData[1024])
{
  switch(partitionInfo->dataType)
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


nodeptr pickRandomSubtree(pllInstance *tr)
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

  
void computeAllAncestralVectors(nodeptr p, pllInstance *tr, partitionList *pr)
{
  /* if this is not a tip, for which evidently it does not make sense 
     to compute the ancestral sequence because we have the real one ....
  */

  if(!isTip(p->number, tr->mxtips))
    {
      /* descend recursively to compute the ancestral states in the left and right subtrees */

      computeAllAncestralVectors(p->next->back, tr, pr);
      computeAllAncestralVectors(p->next->next->back, tr, pr);
      
      /* then compute the ancestral state at node p */

      newviewGenericAncestral(tr, pr, p);

      /* and print it to terminal, the two booleans that are set to PLL_TRUE here 
	 tell the function to print the marginal probabilities as well as 
	 a discrete inner sequence, that is, ACGT etc., always selecting and printing 
	 the state that has the highest probability */

      printAncestralState(p, PLL_TRUE, PLL_TRUE, tr, pr);
    }
}



void initializePartitionData(pllInstance *localTree, partitionList * localPartitions)
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

  int 
    tid  = localTree->threadID,
    innerNodes = localTree->mxtips - 2;

  if(tid > 0)
      localTree->rateCategory    = (int *)    rax_calloc((size_t)localTree->originalCrunchedLength, sizeof(int));	    

  for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
    {
      size_t 
	width = localPartitions->partitionData[model]->width;

      const partitionLengths 
	*pl = getPartitionLengths(localPartitions->partitionData[model]);

      /* 
	 globalScaler needs to be 2 * localTree->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
      */

      localPartitions->partitionData[model]->globalScaler    = (unsigned int *)rax_calloc(2 *(size_t)localTree->mxtips, sizeof(unsigned int));

      localPartitions->partitionData[model]->left              = (double *)rax_malloc_aligned((size_t)pl->leftLength * (maxCategories + 1) * sizeof(double));
      localPartitions->partitionData[model]->right             = (double *)rax_malloc_aligned((size_t)pl->rightLength * (maxCategories + 1) * sizeof(double));
      localPartitions->partitionData[model]->EIGN              = (double*)rax_malloc((size_t)pl->eignLength * sizeof(double));
      localPartitions->partitionData[model]->EV                = (double*)rax_malloc_aligned((size_t)pl->evLength * sizeof(double));
      localPartitions->partitionData[model]->EI                = (double*)rax_malloc((size_t)pl->eiLength * sizeof(double));

      localPartitions->partitionData[model]->substRates        = (double *)rax_malloc((size_t)pl->substRatesLength * sizeof(double));
      localPartitions->partitionData[model]->frequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->freqExponents     = (double*)malloc(pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->empiricalFrequencies       = (double*)rax_malloc((size_t)pl->frequenciesLength * sizeof(double));
      localPartitions->partitionData[model]->tipVector         = (double *)rax_malloc_aligned((size_t)pl->tipVectorLength * sizeof(double));
      
       if(localPartitions->partitionData[model]->dataType == AA_DATA && localPartitions->partitionData[model]->protModels == LG4)      
	{	  	  
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {	    
	      localPartitions->partitionData[model]->EIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));
	      localPartitions->partitionData[model]->EV_LG4[k]                = (double*)rax_malloc_aligned(pl->evLength * sizeof(double));
	      localPartitions->partitionData[model]->EI_LG4[k]                = (double*)rax_malloc(pl->eiLength * sizeof(double));
	      localPartitions->partitionData[model]->substRates_LG4[k]        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));
	      localPartitions->partitionData[model]->frequencies_LG4[k]       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
	      localPartitions->partitionData[model]->tipVector_LG4[k]         = (double *)rax_malloc_aligned(pl->tipVectorLength * sizeof(double));
	    }
	}

      localPartitions->partitionData[model]->symmetryVector    = (int *)rax_malloc((size_t)pl->symmetryVectorLength  * sizeof(int));
      localPartitions->partitionData[model]->frequencyGrouping = (int *)rax_malloc((size_t)pl->frequencyGroupingLength  * sizeof(int));

      localPartitions->partitionData[model]->perSiteRates      = (double *)rax_malloc(sizeof(double) * maxCategories);

      localPartitions->partitionData[model]->nonGTR = PLL_FALSE;

      localPartitions->partitionData[model]->gammaRates = (double*)rax_malloc(sizeof(double) * 4);
      localPartitions->partitionData[model]->yVector = (unsigned char **)rax_malloc(sizeof(unsigned char*) * ((size_t)localTree->mxtips + 1));


      localPartitions->partitionData[model]->xVector = (double **)rax_calloc(sizeof(double*), (size_t)localTree->mxtips);


      /* 
	 Initializing the xVector array like this is absolutely required !!!!
	 I don't know which programming genious removed this, but it must absolutely stay in here!!!!
      */
      
      {
	int k;
	
	for(k = 0; k < localTree->mxtips; k++)
	      localPartitions->partitionData[model]->xVector[k] = (double*)NULL;       
      }


      localPartitions->partitionData[model]->xSpaceVector = (size_t *)rax_calloc((size_t)localTree->mxtips, sizeof(size_t));

      localPartitions->partitionData[model]->sumBuffer = (double *)rax_malloc_aligned(width *
										      (size_t)(localPartitions->partitionData[model]->states) *
										      discreteRateCategories(localTree->rateHetModel) *
										      sizeof(double));

      /* Initialize buffers to store per-site log likelihoods */

      localPartitions->partitionData[model]->perSiteLikelihoods = (double *)rax_malloc_aligned(  width * sizeof(double));

      /* initialize data structures for per-site likelihood scaling */

      if(localTree->fastScaling)
	{
	   localPartitions->partitionData[model]->expVector      = (int **)NULL;
	   localPartitions->partitionData[model]->expSpaceVector = (size_t *)NULL;
	}
      else
	{	 
	  localPartitions->partitionData[model]->expVector      = (int **)rax_malloc(sizeof(int*) * innerNodes);
	   
	  /* 
	     Initializing the expVector array like this is absolutely required !!!!
	     Not doing this can (and did) cause segmentation faults !!!!
	  */
	  
	  {
	    int k;

	    for(k = 0; k < innerNodes; k++)
	      localPartitions->partitionData[model]->expVector[k] = (int*)NULL; 
	  }

	  localPartitions->partitionData[model]->expSpaceVector = (size_t *)rax_calloc(innerNodes, sizeof(size_t));
	}

      /* data structure to store the marginal ancestral probabilities in the sequential version or for each thread */

      localPartitions->partitionData[model]->ancestralBuffer = (double *)rax_malloc_aligned(width *
										 (size_t)(localPartitions->partitionData[model]->states) *
										 sizeof(double));

      /* count and accumulate how many bytes we will need for storing a full ancestral vector. for this we addf over the per-partition space requirements in bytes */
      /* ancestralVectorWidth += ((size_t)(pr->partitionData[model]->upper - pr->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double)); */
      ancestralVectorWidth += ((size_t)(localPartitions->partitionData[model]->upper - localPartitions->partitionData[model]->lower) * (size_t)(localPartitions->partitionData[model]->states) * sizeof(double));
      /* :TODO: do we have to use the original tree for that   */

      localPartitions->partitionData[model]->wgt = (int *)rax_malloc_aligned(width * sizeof(int));

      /* rateCategory must be assigned using rax_calloc() at start up there is only one rate category 0 for all sites */

      localPartitions->partitionData[model]->rateCategory = (int *)rax_calloc(width, sizeof(int));

      if(width > 0 && localTree->saveMemory)
	{
	  localPartitions->partitionData[model]->gapVectorLength = ((int)width / 32) + 1;
	  assert(4 == sizeof(unsigned int));
	  localPartitions->partitionData[model]->gapVector = (unsigned int*)rax_calloc((size_t)localPartitions->partitionData[model]->gapVectorLength * 2 * (size_t)localTree->mxtips, sizeof(unsigned int));
	  localPartitions->partitionData[model]->gapColumn = (double *)rax_malloc_aligned(((size_t)localTree->mxtips) *
									       ((size_t)(localPartitions->partitionData[model]->states)) *
									       discreteRateCategories(localTree->rateHetModel) * sizeof(double));
	}
      else
	{
	  localPartitions->partitionData[model]->gapVectorLength = 0;
	  localPartitions->partitionData[model]->gapVector = (unsigned int*)NULL;
	  localPartitions->partitionData[model]->gapColumn = (double*)NULL;
	}              
    }
}

int virtual_width( int n ) {
    const int global_vw = 2;
    return (n+1) / global_vw * global_vw;
}


void initMemorySavingAndRecom(pllInstance *tr, partitionList *pr)
{
  pllInstance  
    *localTree = tr; 
  partitionList
    *localPartitions = pr;
  size_t model; 

  /* initialize gap bit vectors at tips when memory saving option is enabled */

  if(localTree->saveMemory)
    {
      for(model = 0; model < (size_t)localPartitions->numberOfPartitions; model++)
	{
	  int        
	    undetermined = getUndetermined(localPartitions->partitionData[model]->dataType);

	  size_t
	    i,
	    j,
	    width =  localPartitions->partitionData[model]->width;

	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(localTree->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(localPartitions->partitionData[model]->yVector[j][i] == undetermined)
		    localPartitions->partitionData[model]->gapVector[localPartitions->partitionData[model]->gapVectorLength * j + i / 32] |= mask32[i % 32];
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

double get_branch_length(pllInstance *tr, nodeptr p, int partition_id)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z = p->z[partition_id];
  if(z < PLL_ZMIN) z = PLL_ZMIN;
  if(z > PLL_ZMAX) z = PLL_ZMAX;
  return (-log(z) * tr->fracchange);
}
void set_branch_length(pllInstance *tr, nodeptr p, int partition_id, double bl)
{
  //assert(partition_id < tr->numBranches);
  assert(partition_id < NUM_BRANCHES);
  assert(partition_id >= 0);
  assert(tr->fracchange != -1.0);
  double z;
  z = exp((-1 * bl)/tr->fracchange);
  if(z < PLL_ZMIN) z = PLL_ZMIN;
  if(z > PLL_ZMAX) z = PLL_ZMAX;
  p->z[partition_id] = z;
}

void initializePartitionsSequential(pllInstance *tr, partitionList *pr)
{ 
  size_t
    model;

  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
    assert(pr->partitionData[model]->width == pr->partitionData[model]->upper - pr->partitionData[model]->lower);

  initializePartitionData(tr, pr);

  /* figure in tip sequence data per-site pattern weights */ 
  for(model = 0; model < (size_t)pr->numberOfPartitions; model++)
  {
    size_t
      j;
    size_t lower = pr->partitionData[model]->lower;
    size_t width = pr->partitionData[model]->upper - lower;

    for(j = 1; j <= (size_t)tr->mxtips; j++)
    {
      pr->partitionData[model]->yVector[j] = &(tr->yVector[j][pr->partitionData[model]->lower]);
    }

    memcpy((void*)(&(pr->partitionData[model]->wgt[0])),         (void*)(&(tr->aliaswgt[lower])),      sizeof(int) * width);
  }  

  initMemorySavingAndRecom(tr, pr);
}


/* interface to outside  */
void initializePartitions(pllInstance *tr, pllInstance *localTree, partitionList *pr, partitionList *localPr, int tid, int n)
{
#if (defined(_FINE_GRAIN_MPI) || defined(_USE_PTHREADS))
  initializePartitionsMaster(tr,localTree,pr,localPr,tid,n);
#else
  initializePartitionsSequential(tr, pr);
#endif
}

void 
pllPartitionsDestroy (partitionList ** partitions, int models, int tips)
{
  int i, j;
  partitionList * pl = *partitions;

  for (i = 0; i < models; ++ i)
   {
     rax_free (pl->partitionData[i]->gammaRates);
     rax_free (pl->partitionData[i]->perSiteRates);
     rax_free (pl->partitionData[i]->globalScaler);
     rax_free (pl->partitionData[i]->left);
     rax_free (pl->partitionData[i]->right);
     rax_free (pl->partitionData[i]->EIGN);
     rax_free (pl->partitionData[i]->EV);
     rax_free (pl->partitionData[i]->EI);
     rax_free (pl->partitionData[i]->substRates);
     rax_free (pl->partitionData[i]->frequencies);
     rax_free (pl->partitionData[i]->empiricalFrequencies);
     rax_free (pl->partitionData[i]->tipVector);
     rax_free (pl->partitionData[i]->symmetryVector);
     rax_free (pl->partitionData[i]->frequencyGrouping);
     for (j = 0; j < tips; ++ j)
       rax_free (pl->partitionData[i]->xVector[j]);
     rax_free (pl->partitionData[i]->xVector);
     rax_free (pl->partitionData[i]->yVector);
     rax_free (pl->partitionData[i]->xSpaceVector);
     rax_free (pl->partitionData[i]->sumBuffer);
     rax_free (pl->partitionData[i]->ancestralBuffer);
     rax_free (pl->partitionData[i]->wgt);
     rax_free (pl->partitionData[i]->rateCategory);
     rax_free (pl->partitionData[i]->gapVector);
     rax_free (pl->partitionData[i]->gapColumn);
     rax_free (pl->partitionData[i]->perSiteLikelihoods);
     rax_free (pl->partitionData[i]->partitionName);
     rax_free (pl->partitionData[i]->expSpaceVector);
     /*TODO: Deallocate all entries of expVector */
     if (pl->partitionData[i]->expVector)
      {
        for (j = 0; j < tips - 2; ++ j)
          rax_free (pl->partitionData[i]->expVector[j]);
      }
     rax_free (pl->partitionData[i]->expVector);
     rax_free (pl->partitionData[i]);
   }
  rax_free (pl->partitionData);
  rax_free (pl);

  *partitions = NULL;
}

/** @brief Correspondance check between partitions and alignment

    This function checks whether the partitions to be created and the given
    alignment correspond, that is, whether each site of the alignment is
    assigned to exactly one partition.

    @param parts
      A list of partitions suggested by the caller

    @param phylip
      The multiple sequence alignment
    
    @return
      Returns \a 1 in case of success, otherwise \a 0
*/
int
pllPartitionsValidate (struct pllQueue * parts, struct pllPhylip * phylip)
{
  int nparts;
  char * used;
  struct pllQueueItem * elm;
  struct pllQueueItem * regionItem;
  struct pllPartitionRegion * region;
  struct pllPartitionInfo * pi;
  int i;

  /* check if the list contains at least one partition */
  nparts = pllQueueSize (parts);
  if (!nparts) return (0);


  /* boolean array for marking that a site was assigned a partition */
  used = (char *) rax_calloc (phylip->seqLen, sizeof (char));

  /* traverse all partitions and their respective regions and mark sites */
  for (elm = parts->head; elm; elm = elm->next)
   {
     pi = (struct pllPartitionInfo *) elm->item;
     
     for (regionItem = pi->regionList->head; regionItem; regionItem = regionItem->next)
      {
        region = (struct pllPartitionRegion *) regionItem->item;

        for (i = region->start - 1; i < region->end; i += region->stride)
         {
           if (used[i])
            {
              rax_free (used);
              return (0);
            }
           used[i] = 1; 
         }
      }
   }

  /* check whether all sites were assigned a partition */
  for (i = 0; i < phylip->seqLen; ++ i)
    if (used[i] != 1)
     {
       rax_free (used);
       return (0);
     }

  rax_free (used);
  return (1);
}

/** @brief Swap two sites in a buffer
    
    Swaps sites \a s1 and \a s2 in buffer \a buf which consists of \a nTaxa + 1
    taxa (i.e. rows), and the first row contains no information, i.e. it is not
    accessed.

    @param buffer
      Memory buffer

    @param s1
      First site

    @param s2
      Second site

    @param nTaxa
      Number of taxa, i.e. size of site
*/
static inline void
swapSite (unsigned char ** buf, int s1, int s2, int nTaxa)
{
  int i;
  int x;

  for (i = 1; i <= nTaxa; ++ i)
  {
    x = buf[i][s1];
    buf[i][s1] = buf[i][s2];
    buf[i][s2] = x;
  }
}

/** @brief Constructs the list of partitions according to the proposed partition scheme
    
    A static function that construcs the \a partitionList structure according to
    the partition scheme \b AFTER the sites have been repositioned in contiguous
    regions according to the partition scheme.

    @param bounds
      An array of the new starting and ending posititons of sites in the alignment for each partition.
      This array is of size 2 * \a nparts. The elements are always couples (lower,upper). The upper
      bounds is a site that is not included in the partition

    @todo
      Fix the bug in PLL 

    @param nparts
      The number of partitions to be created
      
*/
static partitionList *
createPartitions (struct pllQueue * parts, int * bounds)
{
  partitionList * pl;
  struct pllPartitionInfo * pi;
  struct pllQueueItem * elm;
  int i;

  pl = (partitionList *) rax_malloc (sizeof (partitionList));
  
  // TODO: fix this
  pl->perGeneBranchLengths =      0;

  // TODO: change NUM_BRANCHES to number of partitions I guess
  pl->partitionData = (pInfo **) rax_calloc (NUM_BRANCHES, sizeof (pInfo *));
  
  for (i = 0, elm = parts->head; elm; elm = elm->next, ++ i)
   {
     pi = (struct pllPartitionInfo *) elm->item;
     pl->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));

     pl->partitionData[i]->lower = bounds[i << 1];
     pl->partitionData[i]->upper = bounds[(i << 1) + 1];
     pl->partitionData[i]->width = bounds[(i << 1) + 1] - bounds[i << 1];

     if (pi->dataType == DNA_DATA)
      {
        pl->partitionData[i]->protModels                = -1;
        pl->partitionData[i]->protFreqs                 = -1;
        pl->partitionData[i]->dataType                  = DNA_DATA;
        pl->partitionData[i]->maxTipStates              = 16;
        pl->partitionData[i]->optimizeBaseFrequencies   = pi->optimizeBaseFrequencies;
      }
     else /* AA_DATA */
      {
        pl->partitionData[i]->dataType                  = AA_DATA; 
        pl->partitionData[i]->protModels                = pi->protModels;
        pl->partitionData[i]->maxTipStates              = 23;
        pl->partitionData[i]->protFreqs                 = pi->protFreqs;
        pl->partitionData[i]->protModels                = pi->protModels;
        pl->partitionData[i]->optimizeBaseFrequencies   = pi->optimizeBaseFrequencies;
      }
     
     pl->partitionData[i]->states                = pLengths[pl->partitionData[i]->dataType].states;
     pl->partitionData[i]->numberOfCategories    =        1;
     pl->partitionData[i]->autoProtModels        =        0;
     pl->partitionData[i]->nonGTR                =        0;
     pl->partitionData[i]->protFreqs             =        0;
     pl->partitionData[i]->partitionContribution =     -1.0;
     pl->partitionData[i]->partitionLH           =      0.0;
     pl->partitionData[i]->fracchange            =      1.0;
     pl->partitionData[i]->executeModel          =     PLL_TRUE;


     pl->partitionData[i]->partitionName         = (char *) rax_malloc ((strlen (pi->partitionName) + 1) * sizeof (char));
     strcpy (pl->partitionData[i]->partitionName, pi->partitionName);
   }

  return (pl);
}


/** @brief Constructs the proposed partition scheme 

    This function constructs the proposed partition scheme. It assumes
    that the partition scheme is correct.

    @note This function \b does \b not validate the partition scheme.
    The user must manually call the \fn pllPartitionsValidate function
    for validation
    
    @param parts
      A list of partitions suggested by the caller

    @param phylip
      The multiple sequence alignment

    @return
      Returns a pointer to \a partitionList structure of partitions in case of success, \b NULL otherwise
*/
partitionList *
pllPartitionsCommit (struct pllQueue * parts, struct pllPhylip * phylip)
{
  int * oi;
  int i, j, dst;
  struct pllQueueItem * elm;
  struct pllQueueItem * regionItem;
  struct pllPartitionRegion * region;
  struct pllPartitionInfo * pi;
  partitionList * pl;
  int * newBounds;
  int k, nparts;

  dst = k = 0;
  oi  = (int *) rax_malloc (phylip->seqLen * sizeof (int));
  for (i = 0; i < phylip->seqLen; ++ i) oi[i] = i;

  nparts = pllQueueSize (parts);
  newBounds = (int *) rax_malloc (2 * nparts * sizeof (int));

  /* reposition the sites in the alignment */
  for (elm = parts->head; elm; elm = elm->next, ++ k)
   {
     pi = (struct pllPartitionInfo *) elm->item;
     
     newBounds[k << 1] = dst;   /* set the lower column for this partition */
     for (regionItem = pi->regionList->head; regionItem; regionItem = regionItem->next)
      {
        region = (struct pllPartitionRegion *) regionItem->item;

        for (i = region->start - 1; i < region->end; i += region->stride)
         {
           if (oi[i] == i)
            {
              swapSite (phylip->seq, dst, i, phylip->nTaxa);
              oi[dst++] = i;
            }
           else if (oi[i] < i)
            {
              j = oi[i];
              while (j < i) j = oi[j];

              swapSite (phylip->seq, dst, j, phylip->nTaxa);
              oi[dst++] = j;
            }
         }
      }
     newBounds[(k << 1) + 1] = dst;    /* set the uppwer limit for this partition */
   }
  pl = createPartitions (parts, newBounds);
  pl->numberOfPartitions = nparts;

  rax_free (newBounds);
  rax_free (oi);

  return (pl);
}

/** @brief Copy a site to another buffer

    Copies site \a from from buffer \a src to \a to in buffer \a dst. Both buffers
    must consist of \a nTaxa + 1 taxa and the first row contains no information, i.e.
    it is not accessed.

    @param dst
      Destination buffer

    @param src
      Source buffer

    @param to
      At which position in \a dst to copy the site to

    @param from
      Which site from \a src to copy

    @param nTaxa
      Number of taxa, i.e. size of site
*/
static inline void
copySite (unsigned char ** dst, unsigned char ** src, int to, int from, int nTaxa)
{
  int i;

  for (i = 1; i <= nTaxa; ++ i)
   {
     dst[i][to] = src[i][from];
   }
}

/** @brief Remove duplicate sites from alignment and update weights vector

    Removes duplicate sites from the alignment given the partitions list
    and updates the weight vector of the alignment and the boundaries
    (upper, lower, width) for each partition.

    @param phylip
      The phylip alignment
    
    @param pl
      List of partitions

*/
void 
pllPhylipRemoveDuplicate (struct pllPhylip * phylip, partitionList * pl)
{
  int i, j, k, p;
  char *** sites;
  void ** memptr;
  int ** oi;
  int dups = 0;
  int lower;

  /* allocate space for the transposed alignments (sites) for every partition */
  sites  = (char ***) rax_malloc (pl->numberOfPartitions * sizeof (char **));
  memptr = (void **)  rax_malloc (pl->numberOfPartitions * sizeof (void *));
  oi     = (int **)   rax_malloc (pl->numberOfPartitions * sizeof (int *));

  /* transpose the sites by partition */
  for (p = 0; p < pl->numberOfPartitions; ++ p)
   {
     sites[p]  = (char **) rax_malloc (pl->partitionData[p]->width * sizeof (char *));
     memptr[p] = rax_malloc ((phylip->nTaxa + 1) * pl->partitionData[p]->width * sizeof (char));

     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        sites[p][i] = (char *) (memptr[p] + i * (phylip->nTaxa + 1) * sizeof (char));
      }

     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        for (j = 0; j < phylip->nTaxa; ++ j)
         {
           sites[p][i][j] = phylip->seq[j + 1][pl->partitionData[p]->lower + i]; 
         }
        sites[p][i][j] = 0;
      }

     oi[p] = ssort1main (sites[p], pl->partitionData[p]->width);

     for (i = 0; i < pl->partitionData[p]->width; ++ i) oi[p][i] = 1;

     for (i = 1; i < pl->partitionData[p]->width; ++ i)
      {
        if (! strcmp (sites[p][i], sites[p][i - 1]))
         {
           ++ dups;
           oi[p][i] = 0;
         }
      }
   }

  /* allocate memory for the alignment without duplicates*/
  rax_free (phylip->seq[1]);
  rax_free (phylip->weights);

  phylip->seqLen = phylip->seqLen - dups;
  phylip->seq[0] = (unsigned char *) rax_malloc ((phylip->seqLen + 1) * sizeof (unsigned char) * phylip->nTaxa);
  for (i = 0; i < phylip->nTaxa; ++ i)
   {
     phylip->seq[i + 1] = (unsigned char *) (phylip->seq[0] + i * (phylip->seqLen + 1) * sizeof (unsigned char));
     phylip->seq[i + 1][phylip->seqLen] = 0;
   }

  phylip->weights    = (int *) rax_malloc ((phylip->seqLen) * sizeof (int));
  phylip->weights[0] = 1;

  /* transpose sites back to alignment */
  for (p = 0, k = 0; p < pl->numberOfPartitions; ++ p)
   {
     lower = k;
     for (i = 0; i < pl->partitionData[p]->width; ++ i)
      {
        if (!oi[p][i])
         {
           ++ phylip->weights[k - 1];
         }
        else
         {
           phylip->weights[k] = 1;
           for (j = 0; j < phylip->nTaxa; ++ j)
            {
              phylip->seq[j + 1][k] = sites[p][i][j];
            }
           ++ k;
         }
      }
     pl->partitionData[p]->lower = lower;
     pl->partitionData[p]->upper = k;
     pl->partitionData[p]->width = k - lower;
   }

  /* deallocate storage for transposed alignment (sites) */
  for (p = 0; p < pl->numberOfPartitions; ++ p)
   {
     rax_free (oi[p]);
     rax_free (memptr[p]);
     rax_free (sites[p]);
   }
  rax_free (oi);
  rax_free (sites);
  rax_free (memptr);
}


static void
genericBaseFrequencies (const int numFreqs, struct pllPhylip * phylip, int lower, int upper, boolean smoothFrequencies, const unsigned int * bitMask, double * pfreqs)
{
  double 
    wj, 
    acc,
    sumf[64],   
    temp[64];
 
  int     
    i, 
    j, 
    k, 
    l;

  unsigned char  *yptr;  

  for(l = 0; l < numFreqs; l++)	    
    pfreqs[l] = 1.0 / ((double)numFreqs);
	  
  for (k = 1; k <= 8; k++) 
    {	     	   	    	      			
      for(l = 0; l < numFreqs; l++)
	sumf[l] = 0.0;
	      
      for (i = 1; i <= phylip->nTaxa; i++) 
	{		 
          yptr = phylip->seq[i];
	  
	  for(j = lower; j < upper; j++) 
	    {
	      unsigned int code = bitMask[yptr[j]];
	      assert(code >= 1);
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if((code >> l) & 1)
		    temp[l] = pfreqs[l];
		  else
		    temp[l] = 0.0;
		}		      	      
	      
	      for(l = 0, acc = 0.0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)
		    acc += temp[l];
		}
	      
	      wj = phylip->weights[j] / acc;
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)		    
		    sumf[l] += wj * temp[l];			     				   			     		   
		}
	    }
	}	    	      
      
      for(l = 0, acc = 0.0; l < numFreqs; l++)
	{
	  if(sumf[l] != 0.0)
	    acc += sumf[l];
	}
	      
      for(l = 0; l < numFreqs; l++)
	pfreqs[l] = sumf[l] / acc;	     
    }

   /* TODO: What is that? */
/*
  if(smoothFrequencies)         
   {;
    smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));	   
   }
  else    
    {
      boolean 
	zeroFreq = PLL_FALSE;

      char 
	typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);  

      for(l = 0; l < numFreqs; l++)
	{
	  if(pfreqs[l] == 0.0)
	    {
	      printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
	      printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
	      zeroFreq = PLL_TRUE;
	    }
	}

      if(zeroFreq)
	exit(-1);

      for(l = 0; l < numFreqs; l++)
	{
	  assert(pfreqs[l] > 0.0);
	  tr->partitionData[model].frequencies[l] = pfreqs[l];
	}     
    }  
*/

  
}

double **
pllBaseFrequenciesGTR (partitionList * pl, struct pllPhylip * phylip)
{
  int
    model,
    lower,
    upper,
    states;

  double ** freqs;

  freqs = (double **) rax_malloc (pl->numberOfPartitions * sizeof (double *));

  for (model = 0; model < pl->numberOfPartitions; ++ model)
   {
     lower    = pl->partitionData[model]->lower;
     upper    = pl->partitionData[model]->upper;
     states   = pl->partitionData[model]->states;
     freqs[model] = (double *) rax_malloc (states * sizeof (double));

     switch  (pl->partitionData[model]->dataType)
      {
        case AA_DATA:
        case DNA_DATA:
          genericBaseFrequencies (states, 
                                  phylip, 
                                  lower, 
                                  upper, 
                                  pLengths[pl->partitionData[model]->dataType].smoothFrequencies,
                                  pLengths[pl->partitionData[model]->dataType].bitVector,
                                  freqs[model]
                                 );
          break;
      }
   }

  return (freqs);
}
/*
double ** pllBaseFrequenciesGTR(rawdata *rdta, cruncheddata *cdta, pllInstance *tr)
{  
  int 
    model,
    lower,
    upper,
    states;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      lower = tr->partitionData[model].lower;
      upper = tr->partitionData[model].upper;	  	 
      states = tr->partitionData[model].states;
	
      switch(tr->partitionData[model].dataType)
	{
	case GENERIC_32:
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	    case MK_MULTI_STATE:	   
	      {	       
		int i;
		double 
		  freq = 1.0 / (double)states,
		  acc = 0.0;

		for(i = 0; i < states; i++)
		  {
		    acc += freq;
		    tr->partitionData[model].frequencies[i] = freq;
		    // printf("%f \n", freq);
		  }
		// printf("Frequency Deviation: %1.60f\n", acc);
	      }
	      break;
	     case GTR_MULTI_STATE:
	      genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, PLL_TRUE,
				     bitVector32);
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case GENERIC_64:	 
	  assert(0);
	  break;
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case SECONDARY_DATA:
	case AA_DATA:
	case DNA_DATA:
	case BINARY_DATA:	  
	  genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, 
				 getSmoothFreqs(tr->partitionData[model].dataType),
				 getBitVector(tr->partitionData[model].dataType));	  	 
	  break;	
	default:
	  assert(0);     
	}      
    }
  
  return;
}
*/

void
pllEmpiricalFrequenciesDestroy (double *** empiricalFrequencies, int models)
{
  int i;

  for (i = 0; i < models; ++ i)
   {
     rax_free ((*empiricalFrequencies)[i]);
   }
  rax_free (*empiricalFrequencies);

  *empiricalFrequencies = NULL;
}

/** @brief Load alignment to the PLL instance
    
    Loads the parsed phylip alignment to the PLL instance.
    In case the deep switch is specified, the phylip structure
    will be used as the alignment.

    @param tr
      The library instance

    @param phylip
      The phylip alignment

    @param partitions
      List of partitions

    @param bDeep
      Controls how the alignment is used within the library instance.
      If PLL_DEEP_COPY is specified, new memory will be allocated
      and the alignment will be copied there (deep copy). On the other
      hand, if PLL_SHALLOW_COPY is specified, only the pointers will be
      copied and therefore, the alignment will be shared between the 
      phylip structure and the library instance (shallow copy).

    @return
      Returns 1 in case of success, 0 otherwise.
*/
int
pllLoadAlignment (pllInstance * tr, struct pllPhylip * phylip, partitionList * partitions, int bDeep)
{
  int i;
  nodeptr node;

  if (tr->mxtips != phylip->nTaxa) return (0);

  /* Do the base substitution (from A,C,G....  ->   0,1,2,3....)*/
  pllBaseSubstitute (phylip, partitions);

  tr->aliaswgt = (int *) rax_malloc (phylip->seqLen * sizeof (int));
  memcpy (tr->aliaswgt, phylip->weights, phylip->seqLen * sizeof (int));

  tr->originalCrunchedLength = phylip->seqLen;
  tr->rateCategory           = (int *)   rax_calloc (tr->originalCrunchedLength, sizeof (int));
  tr->patrat                 = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->patratStored           = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));
  tr->lhs                    = (double*) rax_malloc((size_t)tr->originalCrunchedLength * sizeof(double));

  /* allocate memory for the alignment */
  tr->yVector    = (unsigned char **) rax_malloc ((phylip->nTaxa + 1) * sizeof (unsigned char *));                                                                                                                                                                      
  tr->bDeep = bDeep;

  if (bDeep == PLL_DEEP_COPY)
   {
     tr->yVector[0] = (unsigned char *)  rax_malloc (sizeof (unsigned char) * (phylip->seqLen + 1) * phylip->nTaxa);
     for (i = 1; i <= phylip->nTaxa; ++ i)                      
      {                     
        tr->yVector[i] = (unsigned char *) (tr->yVector[0] + (i - 1) * (phylip->seqLen + 1) * sizeof (unsigned char));
        tr->yVector[i][phylip->seqLen] = 0;
      }                     
   }
                         
  /* place sequences to tips */                              
  for (i = 1; i <= phylip->nTaxa; ++ i)                      
   {                     
     if (!pllHashSearch (tr->nameHash, phylip->label[i],(void **)&node)) 
      {
        //rax_free (tr->originalCrunchedLength);
        rax_free (tr->rateCategory);
        rax_free (tr->patrat);
        rax_free (tr->patratStored);
        rax_free (tr->lhs);
        if (bDeep == PLL_DEEP_COPY) rax_free (tr->yVector[0]);
        rax_free (tr->yVector);
        return (0);
      }
     if (bDeep == PLL_DEEP_COPY)
       memcpy (tr->yVector[node->number], phylip->seq[i], phylip->seqLen );
     else
       tr->yVector[node->number] = phylip->seq[i];
   }

  return (1);
}

/** @brief Create the main instance of PLL
    
    Create an instance of the phylogenetic likelihood library

    @param rateHetModel
      Rate heterogeneity model

    @param fastScaling
      explain fastScaling here

    @param saveMemory
      explain saveMemory here

    @param useRecom
      explain useRecom
    
    @todo
      Document fastScaling, rate heterogeneity and saveMemory and useRecom

    @return
      On success returns an instance to PLL, otherwise \b NULL
*/
pllInstance *
pllCreateInstance (int rateHetModel, int fastScaling, int saveMemory, int useRecom, long randomNumberSeed)
{
  pllInstance * tr;

  if (rateHetModel != GAMMA && rateHetModel != CAT) return NULL;

  tr = (pllInstance *) rax_calloc (1, sizeof (pllInstance));

  tr->threadID     = 0;
  tr->rateHetModel = rateHetModel;
  tr->fastScaling  = fastScaling;
  tr->saveMemory   = saveMemory;
  tr->useRecom     = useRecom;
  
  tr->randomNumberSeed = randomNumberSeed;

  /* remove it from the library */
  tr->useMedian    = PLL_FALSE;

  tr->maxCategories = (rateHetModel == GAMMA) ? 4 : 25;
  
  return (tr);
}

/** @brief Initialize PLL tree structure with default values
    
    Initialize PLL tree structure with default values and allocate 
    memory for its elements.

    @todo
      STILL NOT FINISHED
*/
void pllTreeInitDefaults (pllInstance * tr, int nodes, int tips)
{
  nodeptr p0, p, q;
  int i, j;
  int inner;

  

  /* TODO: make a proper static setupTree function */

  inner = tips - 1;

  tr->mxtips = tips;

  tr->bigCutoff = PLL_FALSE;
  tr->treeStringLength = tr->mxtips * (PLL_NMLNGTH + 128) + 256 + tr->mxtips * 2;
  tr->tree_string = (char *) rax_calloc ( tr->treeStringLength, sizeof(char));
  tr->tree0 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));

  
  p0 = (nodeptr) rax_malloc ((tips + 3 * inner) * sizeof (node));
  assert (p0);

  tr->nodeBaseAddress  = p0;

  tr->nameList         = (char **)   rax_malloc ((tips + 1) * sizeof (char *));
  tr->nodep            = (nodeptr *) rax_malloc ((2 * tips) * sizeof (nodeptr));
  assert (tr->nameList && tr->nodep);

  tr->nodep[0] = NULL;          


  /* TODO: FIX THIS! */
  tr->fracchange = -1;

  for (i = 1; i <= tips; ++ i)
   {
     p = p0++;

     //p->hash      = KISS32();     
     p->x         = 0;
     p->xBips     = 0;
     p->number    = i;
     p->next      = p;
     p->back      = NULL;
     p->bInf      = NULL;
     tr->nodep[i]  = p;
   }

  for (i = tips + 1; i <= tips + inner; ++i)
   {
     q = NULL;
     for (j = 1; j <= 3; ++ j)
     {
       p = p0++;
       if (j == 1)
        {
          p->xBips = 1;
          p->x = 1; //p->x     = 1;
        }
       else
        {
          p->xBips = 0;
          p->x     = 0;
        }
       p->number = i;
       p->next   = q;
       p->bInf   = NULL;
       p->back   = NULL;
       p->hash   = 0;
       q         = p;
     }
    p->next->next->next = p;
    tr->nodep[i]         = p;
   }

  tr->likelihood  = PLL_UNLIKELY;
  tr->start       = NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;

  for (i = 0; i < NUM_BRANCHES; ++ i) tr->partitionSmoothed[i] = PLL_FALSE;

  tr->bitVectors = NULL;
  tr->vLength    = 0;
  tr->h          = NULL;

  /* TODO: Fix hash type */
  tr->nameHash   = pllHashInit (10 * tr->mxtips);

  /* TODO: do these options really fit here or should they be put elsewhere? */
  tr->td[0].count            = 0;
  tr->td[0].ti               = (traversalInfo *) rax_malloc (sizeof(traversalInfo) * (size_t)tr->mxtips);
  tr->td[0].parameterValues  = (double *) rax_malloc(sizeof(double) * (size_t)NUM_BRANCHES);
  tr->td[0].executeModel     = (boolean *) rax_malloc (sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr->td[0].executeModel[0]  = PLL_TRUE;                                                                                                                                                                                                                                    
  for (i = 0; i < NUM_BRANCHES; ++ i) tr->td[0].executeModel[i] = PLL_TRUE;


  
}


/** @brief Initializes the PLL tree topology according to a parsed newick tree

    Set the tree topology based on a parsed and validated newick tree

    @param tree
      The PLL instance

    @param nt
      The \a pllNewickTree wrapper structure that contains the parsed newick tree

    @param bDefaultZ
      If set to \b PLL_TRUE then the branch lengths will be reset to the default
      value.
*/
void
pllTreeInitTopologyNewick (pllInstance * tr, struct pllNewickTree * nt, int bUseDefaultZ)
{
  struct pllStack * nodeStack = NULL;
  struct pllStack * head;
  struct item_t * item;
  int i, j, k;
  
/*
  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     partitions->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));
     partitions->partitionData[i]->partitionContribution = -1.0;
     partitions->partitionData[i]->partitionLH           =  0.0;
     partitions->partitionData[i]->fracchange            =  1.0;
   }
*/
  
  pllTreeInitDefaults (tr, nt->nodes, nt->tips);

  i = nt->tips + 1;
  j = 1;
  nodeptr v;
  
  
  for (head = nt->tree; head; head = head->next)
  {
    item = (struct item_t *) head->item;
    if (!nodeStack)
     {
       pllStackPush (&nodeStack, tr->nodep[i]);
       pllStackPush (&nodeStack, tr->nodep[i]->next);
       pllStackPush (&nodeStack, tr->nodep[i]->next->next);
       ++i;
     }
    else
     {
       v = (nodeptr) pllStackPop (&nodeStack);
       if (item->rank)  /* internal node */
        {
          v->back           = tr->nodep[i];
          tr->nodep[i]->back = v; //t->nodep[v->number]
          pllStackPush (&nodeStack, tr->nodep[i]->next);
          pllStackPush (&nodeStack, tr->nodep[i]->next->next);
          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < NUM_BRANCHES; ++ k)
             v->z[k] = tr->nodep[i]->z[k] = z;

          ++ i;
        }
       else             /* leaf */
        {
          v->back           = tr->nodep[j];
          tr->nodep[j]->back = v; //t->nodep[v->number];

          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < NUM_BRANCHES; ++ k)
            v->z[k] = tr->nodep[j]->z[k] = z;
            
          //t->nameList[j] = strdup (item->name);
          tr->nameList[j] = (char *) rax_malloc ((strlen (item->name) + 1) * sizeof (char));
          strcpy (tr->nameList[j], item->name);
          
          pllHashAdd (tr->nameHash, tr->nameList[j], (void *) (tr->nodep[j]));
          ++ j;
        }
     }
  }
  
  tr->start = tr->nodep[1];
  
  pllStackClear (&nodeStack);

  if (bUseDefaultZ == PLL_TRUE) resetBranches (tr);
}

/** @brief Initialize PLL tree with a random topology

    Initializes the PLL tree with a randomly created topology

    @todo
      Perhaps pass a seed?

    @param tr
      The PLL instance

    @param tips
      Number of tips

    @param nameList
      A set of \a tips names representing the taxa labels
*/
void 
pllTreeInitTopologyRandom (pllInstance * tr, int tips, char ** nameList)
{
  int i;
  pllTreeInitDefaults (tr, 2 * tips - 1, tips);

  for (i = 1; i <= tips; ++ i)
   {
     tr->nameList[i] = (char *) rax_malloc ((strlen (nameList[i]) + 1) * sizeof (char));
     strcpy (tr->nameList[i], nameList[i]);
     pllHashAdd (tr->nameHash, tr->nameList[i], (void *) (tr->nodep[i]));
   }

  makeRandomTree (tr);
}

void
pllBaseSubstitute (struct pllPhylip * phylip, partitionList * partitions)
{
  char meaningDNA[256];
  char  meaningAA[256];
  char * d;
  int i, j, k;

  for (i = 0; i < 256; ++ i)
   {
     meaningDNA[i] = -1;
     meaningAA[i]  = -1;
   }

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['Y'] = 10;
  meaningDNA['a'] =  1;
  meaningDNA['b'] = 14;
  meaningDNA['c'] =  2;
  meaningDNA['d'] = 13;
  meaningDNA['g'] =  4;
  meaningDNA['h'] = 11;
  meaningDNA['k'] = 12;
  meaningDNA['m'] =  3;
  meaningDNA['r'] =  5;
  meaningDNA['s'] =  6;
  meaningDNA['t'] =  8;
  meaningDNA['u'] =  8;
  meaningDNA['v'] =  7;
  meaningDNA['w'] =  9;
  meaningDNA['y'] = 10;

  meaningDNA['N'] =
  meaningDNA['n'] =
  meaningDNA['O'] =
  meaningDNA['o'] =
  meaningDNA['X'] =
  meaningDNA['x'] =
  meaningDNA['-'] =
  meaningDNA['?'] = 15;
 
  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/
  meaningAA['a'] =  0;  /* alanine */
  meaningAA['r'] =  1;  /* arginine */
  meaningAA['n'] =  2;  /*  asparagine*/
  meaningAA['d'] =  3;  /* aspartic */
  meaningAA['c'] =  4;  /* cysteine */
  meaningAA['q'] =  5;  /* glutamine */
  meaningAA['e'] =  6;  /* glutamic */
  meaningAA['g'] =  7;  /* glycine */
  meaningAA['h'] =  8;  /* histidine */
  meaningAA['i'] =  9;  /* isoleucine */
  meaningAA['l'] =  10; /* leucine */
  meaningAA['k'] =  11; /* lysine */
  meaningAA['m'] =  12; /* methionine */
  meaningAA['f'] =  13; /* phenylalanine */
  meaningAA['p'] =  14; /* proline */
  meaningAA['s'] =  15; /* serine */
  meaningAA['t'] =  16; /* threonine */
  meaningAA['w'] =  17; /* tryptophan */
  meaningAA['y'] =  18; /* tyrosine */
  meaningAA['v'] =  19; /* valine */
  meaningAA['b'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
  meaningAA['x'] = 
  meaningAA['?'] = 
  meaningAA['*'] = 
  meaningAA['-'] = 22;

  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     d = (partitions->partitionData[i]->dataType == DNA_DATA) ? meaningDNA : meaningAA;
     
     for (j = 1; j <= phylip->nTaxa; ++ j)
      {
        for (k = partitions->partitionData[i]->lower; k < partitions->partitionData[i]->upper; ++ k)
         {
           phylip->seq[j][k] = d[phylip->seq[j][k]];
         }
      }
   }
}

/** @brief Deallocate the PLL instance

    Deallocates the library instance and all its elements.

    @param tr
      The PLL instance
*/
void
pllTreeDestroy (pllInstance * tr)
{
  int i;
  for (i = 1; i <= tr->mxtips; ++ i)
    rax_free (tr->nameList[i]);
  
  pllHashDestroy (&(tr->nameHash), PLL_FALSE);
  if (tr->yVector)
   {
     if (tr->bDeep == PLL_DEEP_COPY)
      {
        if (tr->yVector[0]) rax_free (tr->yVector[0]);
      }
     rax_free (tr->yVector);
   }
  rax_free (tr->aliaswgt);
  rax_free (tr->rateCategory);
  rax_free (tr->patrat);
  rax_free (tr->patratStored);
  rax_free (tr->lhs);
  rax_free (tr->td[0].parameterValues);
  rax_free (tr->td[0].executeModel);
  rax_free (tr->td[0].ti);
  rax_free (tr->nameList);
  rax_free (tr->nodep);
  rax_free (tr->nodeBaseAddress);
  rax_free (tr->tree_string);
  rax_free (tr->tree0);
  rax_free (tr->tree1);
  rax_free (tr);
}

/** @brief Initialize partitions according to model parameters
    
    Initializes partitions according to model parameters.

    @param tr
      The PLL instance

    @param bEmpiricalFreqs
      Use empirical frequencies
    
    @param bResetBranches
      Reset branch lengths to default lengths

    @param phylip
      The alignment

    @param partitions
      List of partitions

    @todo
      implement the bEmpiricalFreqs flag
*/
void pllInitModel (pllInstance * tr, int bEmpiricalFreqs, struct pllPhylip * phylip, partitionList * partitions)
{
  double ** ef;

  ef = pllBaseFrequenciesGTR (partitions, phylip);
  initializePartitions (tr, tr, partitions, partitions, 0, 0);
  initModel (tr, ef, partitions);
  pllEmpiricalFrequenciesDestroy (&ef, partitions->numberOfPartitions);
}
