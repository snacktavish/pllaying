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
	  printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
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
	  printf("The file %s RAxML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
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

static boolean setupTree (tree *tr)
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
#endif

}




static boolean modelExists(char *model, tree *tr)
{
  /********** BINARY ********************/

   if(strcmp(model, "PSR\0") == 0)
    {
      tr->rateHetModel = CAT;
      return TRUE;
    }

  if(strcmp(model, "GAMMA\0") == 0)
    {
      tr->rateHetModel = GAMMA;
      return TRUE;
    }

  
  return FALSE;
}



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
      printf(": illegal option -- %c \n", c);
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
	      printf(": option requires an argument -- %c\n", c);
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




/*********************************** *********************************************************/


static void printVersionInfo(void)
{
  printf("\n\nThis is %s version %s released by Alexandros Stamatakis on %s.\n\n",  programName, programVersion, programDate); 
}

static void printMinusFUsage(void)
{
  printf("\n");
 

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");
  
  printf("\n");

  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");

  printf("\n");
}


static void printREADME(void)
{
  printVersionInfo();
  printf("\n");  
  printf("\nTo report bugs send an email to raxml@h-its.org\n");
  printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");

  printf("raxmlLight|raxmlLight-PTHREADS|raxmlLight-MPI\n");
  printf("      -s sequenceFileName| -G binarySequnceFile\n");
  printf("      -n outputFileName\n");
  printf("      -m substitutionModel\n");
  printf("      -t userStartingTree| -R binaryCheckpointFile\n");
  printf("      [-B]\n"); 
  printf("      [-c numberOfCategories]\n");
  printf("      [-D]\n");
  printf("      [-e likelihoodEpsilon] \n");
  printf("      [-f d|o]\n");    
  printf("      [-h] \n");
  printf("      [-i initialRearrangementSetting] \n");
  printf("      [-M]\n");
  printf("      [-P proteinModel]\n");
  printf("      [-q multipleModelFileName] \n");
#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))
  printf("      [-Q]\n");
#endif
  printf("      [-S]\n");
  printf("      [-T numberOfThreads]\n");  
  printf("      [-v]\n"); 
  printf("      [-w outputDirectory] \n"); 
  printf("      [-X]\n");
  printf("\n");
  printf("      -B      Parse phylip file and conduct pattern compression, then store the output in a \n");
  printf("              binary file called sequenceFileName.binary that can be read via the \"-G\" option\n");
  printf("              ATTENTION: the \"-B\" option only works with the sequential version\n");
  printf("\n");
  printf("      -c      Specify number of distinct rate catgories for RAxML when modelOfEvolution\n");
  printf("              is set to GTRPSR\n");
  printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
  printf("              categories to accelerate computations. \n");
  printf("\n");
  printf("              DEFAULT: 25\n");
  printf("\n");
  printf("      -D      ML search convergence criterion. This will break off ML searches if the relative \n");
  printf("              Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles\n");
  printf("              is smaller or equal to 1%s. Usage recommended for very large datasets in terms of taxa.\n", "%");
  printf("              On trees with more than 500 taxa this will yield execution time improvements of approximately 50%s\n",  "%");
  printf("              While yielding only slightly worse trees.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");    
  printf("\n");
  printf("      -e      set model optimization precision in log likelihood units for final\n");
  printf("              optimization of model parameters\n");
  printf("\n");
  printf("              DEFAULT: 0.1 \n"); 
  printf("\n");
  printf("      -f      select algorithm:\n");

  printMinusFUsage();
 
  printf("\n");
  printf("      -G      Read in a binary alignment file (instead of a text-based phylip file with \"-s\") that was previsouly\n");
  printf("              generated with the \"-B\" option. This can substantially save time spent in input parsing \n");
  printf("              for very large parallel runs\n");
  printf("\n");
  printf("      -h      Display this help message.\n");
  printf("\n");  
  printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
  printf("              changes phase\n");
  printf("\n");
  printf("      -m      Model of  Nucleotide or Amino Acid Substitution: \n");
  printf("\n"); 
  printf("              NUCLEOTIDES:\n\n");
  printf("                \"-m GTRPSR\"         : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency.\n");
  printf("                \"-m GTRGAMMA\"       : GTR + GAMMA model of rate heterogeneity. This uses 4 hard-coded discrete rates\n");
  printf("                                      to discretize the GAMMA distribution.\n");
  printf("\n");
  printf("              AMINO ACIDS:\n\n");
  printf("                \"-m PROTPSRmatrixName[F]\"         : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.\n");  
  printf("                \"-m PROTGAMMAmatrixName[F]\"       : specified AA matrix + GAMMA model of rate heterogeneity. This uses 4 hard-coded discrete rates\n");
  printf("                                                    to discretize the GAMMA distribution.\n");
  printf("\n");
  printf("                Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA,\n");
  printf("                PMB, HIVB, HIVW, JTTDCMUT, FLU, AUTO, GTR\n");
  printf("                With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
  printf("                Please note that for mixed models you can in addition specify the per-gene AA model in\n");
  printf("                the mixed model file (see manual for details). Also note that if you estimate AA GTR parameters on a partitioned\n");
  printf("                dataset, they will be linked (estimated jointly) across all partitions to avoid over-parametrization\n");
  printf("                When AUTO is used RAxML will conduct an ML estimate of all available pre-defined AA models (excluding GTR) every time the model parameters\n");
  printf("                are optimized during the tree search.\n");
  printf("                WARNING: we have not figured out yet how to best do this for partitioned analyses, so don't use AUTO for datasets with partitions\n");
  printf("\n");
  printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");
  printf("              Branch lengths for individual partitions will be printed to separate files\n");
  printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
  printf("\n"),
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -n      Specifies the name of the output file.\n");
  printf("\n"); 
  printf("      -P      Specify the file name of a user-defined AA (Protein) substitution model. This file must contain\n");
  printf("              420 entries, the first 400 being the AA substitution rates (this must be a symmetric matrix) and the\n");
  printf("              last 20 are the empirical base frequencies\n");
  printf("\n");
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n");
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n");  
  printf("\n");
#if (defined(_USE_PTHREADS) || (_FINE_GRAIN_MPI))
  printf("      -Q      Enable alternative data/load distribution algorithm for datasets with many partitions\n");
  printf("              In particular under PSR this can lead to parallel performance improvements of over 50 per cent\n");
#endif
  printf("\n");
  printf("      -R      read in a binary checkpoint file called RAxML_binaryCheckpoint.RUN_ID_number\n");
  printf("\n");
  printf("      -s      Specify the name of the alignment data file in PHYLIP format\n");
  printf("\n");
  printf("      -S      turn on memory saving option for gappy multi-gene alignments. For large and gappy datasets specify -S to save memory\n");
  printf("              This will produce slightly different likelihood values, may be a bit slower but can reduce memory consumption\n");
  printf("              from 70GB to 19GB on very large and gappy datasets\n");
  printf("\n");
  printf("      -t      Specify a user starting tree file name in Newick format\n");
  printf("\n");
  printf("      -T      PTHREADS VERSION ONLY! Specify the number of threads you want to run.\n");
  printf("              Make sure to set \"-T\" to at most the number of CPUs you have on your machine,\n");
  printf("              otherwise, there will be a huge performance decrease!\n");
  printf("\n");  
  printf("      -v      Display version information\n");
  printf("\n");
  printf("      -w      FULL (!) path to the directory into which RAxML shall write its output files\n");
  printf("\n");
  printf("              DEFAULT: current directory\n"); 
  printf("\n");
  printf("      -X      EXPERIMENTAL OPTION: This option will do a per-site estimate of protein substitution models\n");
  printf("              by looping over all given, fixed models LG, WAG, JTT, etc and using their respective base frequencies to independently\n");
  printf("              assign a prot subst. model to each site via ML optimization\n");
  printf("              At present this option only works with the GTR+GAMMA model, unpartitioned datasets, and in the sequential\n");
  printf("              version only.\n");
  printf("\n\n\n\n");

}




static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    *optarg,
    model[1024] = "",       
    modelChar;

  double 
    likelihoodEpsilon;
  
  int  
    optind = 1,        
    c,
    nameSet = 0,
    treeSet = 0,   
    modelSet = 0;

  boolean 
    byteFileSet = FALSE;


 
 

  /*********** tr inits **************/

  tr->numberOfThreads = 1; 
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->rateHetModel = GAMMA;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;

  tr->manyPartitions = FALSE;
  
  tr->startingTree = randomTree;
  tr->randomNumberSeed = 12345;

  tr->categories             = 25;

  tr->grouped = FALSE;
  tr->constrained = FALSE;

  tr->gapyness               = 0.0; 

  /********* tr inits end*************/




  while(!bad_opt && ((c = mygetopt(argc,argv,"T:R:e:c:f:i:m:r:t:w:n:s:vhMSDQXbp", &optind, &optarg))!=-1))
    {
    switch(c)
      {    
      case 'p':
	tr->startingTree = parsimonyTree;
	break;
      case 'r':
	sscanf(optarg,"%ld", &(tr->randomNumberSeed));	
	if(tr->randomNumberSeed <= 0)
	  {
	    printf("Random number seed specified via -r must be greater than zero\n");
	    exit(-1);
	  }
	break;
      case 'Q':
#ifdef _USE_PTHREADS
	tr->manyPartitions = TRUE;
#else
	printf("The -Q option does not have an effect in the sequential version\n");
	tr->manyPartitions = FALSE;
#endif
	break;
      case 's':		 	
	strcpy(byteFileName, optarg);	 	
	byteFileSet = TRUE;       
	break;      
      case 'S':
	tr->saveMemory = TRUE;
	break;
      case 'D':
	tr->searchConvergenceCriterion = TRUE;	
	break;
      case 'R':
	adef->useCheckpoint = TRUE;
	strcpy(binaryCheckpointInputName, optarg);
	break;          
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;        
      case 'T':
#ifdef _USE_PTHREADS
	sscanf(optarg,"%d", &(tr->numberOfThreads));
#else      
	printf("Option -T does not have any effect with the sequential or parallel MPI version.\n");
	printf("It is used to specify the number of threads for the Pthreads-based parallelization\n");       
#endif
	break;                       
      case 'e':
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;
	break;    
      
      case 'v':
	printVersionInfo();
	exit(0);
      
      case 'h':
	printREADME();
	exit(0);
     
      case 'c':
	sscanf(optarg, "%d", &tr->categories);
	break;     
      case 'f':
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {	 
	  case 'd':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;	  
	  case 'o':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;	    	  	  	     
	  default:
	    {	     	      
	      printf("Error select one of the following algorithms via -f :\n");
	      printMinusFUsage();	      
	      exit(-1);
	    }
	  }
	break;
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;
      case 'w':
        strcpy(resultDir, optarg);
	resultDirSet = TRUE;
        break;
      case 't':
	strcpy(tree_file, optarg);
	tr->startingTree = givenTree;
	treeSet = 1;       
	break;     
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, tr) == 0)
	  {	    
	    printf("Rate heterogeneity Model %s does not exist\n\n", model);               
	    printf("For per site rates (called CAT in previous versions) use: PSR\n");	
	    printf("For GAMMA use: GAMMA\n");			  
	    exit(-1);
	  }
	else
	  modelSet = 1;
	break;
      case 'b':
#ifdef _BAYESIAN 	
	adef->bayesian = TRUE;
	printf("EXPERIMENTAL BAYESIAN ANALYSIS\n");
	break;
#else
	printf("recompile with Bayesian Makefile to use the -b option \n");
	exit(-1);
	break;
#endif
      default:
	exit(-1);
      }
    }



  

  if(!byteFileSet)
    {      
      printf("\n Error, you must specify a binary format data file with the \"-s\" option\n");
      exit(-1);
    }

  if(!modelSet)
    {     
      printf("\n Error, you must specify a model of rate heterogeneity with the \"-m\" option\n");
      exit(-1);
    }

  if(!nameSet)
    {     
      printf("\n Error: please specify a name for this run with -n\n");
      exit(-1);
    }

  
  
  
  
   {
#ifdef WIN32
    const 
      char *separator = "\\";
#else
    const 
      char *separator = "/";
#endif

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	
#ifndef WIN32
	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
#endif
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
   }


  

  return;
}








static void makeFileNames(void)
{
  int infoFileExists = 0;


  
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);  
  strcpy(infoFileName,         workdir);
  strcpy(randomFileName,       workdir);  
  strcpy(binaryCheckpointName, workdir);
   
  strcat(resultFileName,       "RAxML_result.");
  strcat(logFileName,          "RAxML_log.");  
  strcat(infoFileName,         "RAxML_info.");
  strcat(randomFileName,       "RAxML_randomTree.");  
  strcat(binaryCheckpointName, "RAxML_binaryCheckpoint.");
  
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);  
  strcat(infoFileName,         run_id);
  strcat(randomFileName,       run_id);  
  strcat(binaryCheckpointName, run_id);
  



 
  infoFileExists = filexists(infoFileName);
  
  if(infoFileExists)
    {
      printf("RAxML output files with the run ID <%s> already exist \n", run_id);
      printf("in directory %s ...... exiting\n", workdir);
      
      exit(-1);
    }
}




 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  
  int i, model;
  FILE *infoFile = myfopen(infoFileName, "ab");
  char modelType[128];
  
  
  
  strcpy(modelType, "GAMMA");   
  
  printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);
  
  
  
  if(!adef->compressPatterns)
    printBoth(infoFile, "\nAlignment has %d columns\n\n",  tr->originalCrunchedLength);
  else
    printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->originalCrunchedLength);
  
  
  
  printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * tr->gapyness, "%");
  
  
  switch(adef->mode)
    {	
    case  BIG_RAPID_MODE:	 
      printBoth(infoFile, "\nRAxML rapid hill-climbing mode\n\n");
      break;	
    default:
      assert(0);
    }
 
  if(adef->perGeneBranchLengths)
    printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
  else
    printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);	
   
  printBoth(infoFile, "All free model parameters will be estimated by RAxML\n");
      	
  if(tr->rateHetModel == GAMMA)
    printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
  else    
    printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", tr->categories);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      printBoth(infoFile, "Partition: %d\n", model);
      printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
      printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);
      
      switch(tr->partitionData[model].dataType)
	{
	case DNA_DATA:
	  printBoth(infoFile, "DataType: DNA\n");	     
	  printBoth(infoFile, "Substitution Matrix: GTR\n");
	  break;
	case AA_DATA:
	  assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
	  printBoth(infoFile, "DataType: AA\n");	      
	  printBoth(infoFile, "Substitution Matrix: %s\n", tr->partitionData[model].protModels);
	  printBoth(infoFile, "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");	     
	  break;
	case BINARY_DATA:
	  printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");	      
	  printBoth(infoFile, "Substitution Matrix: Uncorrected\n");
	  break;
	case SECONDARY_DATA:
	  printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");	     
	  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	  break;
	case SECONDARY_DATA_6:
	  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");	     
	  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	  break;
	case SECONDARY_DATA_7:
	  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");	      
	  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
	  break;
	case GENERIC_32:
	  printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	      printBoth(infoFile, "Substitution Matrix: Ordered Likelihood\n");
	      break;
	    case MK_MULTI_STATE:
	      printBoth(infoFile, "Substitution Matrix: MK model\n");
	      break;
	    case GTR_MULTI_STATE:
	      printBoth(infoFile, "Substitution Matrix: GTR\n");
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case GENERIC_64:
	  printBoth(infoFile, "DataType: Codon\n");		  
	  break;		
	default:
	  assert(0);
	}
      printBoth(infoFile, "\n\n\n");
    }
  
  printBoth(infoFile, "\n");
  
  printBoth(infoFile, "RAxML was called as follows:\n\n");
  for(i = 0; i < argc; i++)
    printBoth(infoFile,"%s ", argv[i]);
  printBoth(infoFile,"\n\n\n");
  
  fclose(infoFile);
}

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





static void finalizeInfoFile(tree *tr, analdef *adef)
{
  double t;
  
  t = gettime() - masterTime;
  accumulatedTime = accumulatedTime + t;
  
  switch(adef->mode)
    {	
    case  BIG_RAPID_MODE:	 
      printBothOpen("\n\nOverall Time for 1 Inference %f\n", t);
      printBothOpen("\nOverall accumulated Time (in case of restarts): %f\n\n", accumulatedTime);
      printBothOpen("Likelihood   : %f\n", tr->likelihood);
      printBothOpen("\n\n");	  	  
      printBothOpen("Final tree written to:                 %s\n", resultFileName);
      printBothOpen("Execution Log File written to:         %s\n", logFileName);
      printBothOpen("Execution information file written to: %s\n",infoFileName);	
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
     
      memcpy(localTree->td[0].executeModel,    tr->td[0].executeModel,    sizeof(boolean) * localTree->NumberOfModels);
      memcpy(localTree->td[0].parameterValues, tr->td[0].parameterValues, sizeof(double) * localTree->NumberOfModels);
      
      if(localTree->td[0].traversalHasChanged)
	memcpy(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));
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
      makeGammaCats(localTree->partitionData[model].alpha, localTree->partitionData[model].gammaRates, 4); 
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
   
      /* check for errors */
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




  



static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

/*static char bits_in_16bits [0x1u << 16];*/

static void compute_bits_in_16bits(char *bits_in_16bits)
{
    unsigned int i;    
    
    /* size is 65536 */

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);       

    return ;
}

unsigned int precomputed16_bitcount (unsigned int n, char *bits_in_16bits)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}



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
	 Achieving a good balance of alignment sites to cores boils down to the mult-processor scheduling problem known from theoretical comp. sci.
	 which si NP-complete.
	 We have implemented very simple "standard" heuristics for solving the multiprocessor scheduling problem that turn out to work very well
	 and are cheap to compute 
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
    printBothOpen("\nMulti-processor partition data distribution enabled (-Q option)\n");

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
    model,
    maxCategories;

  localTree->threadID = tid; 
  compute_bits_in_16bits(localTree->bits_in_16bits);


#ifdef _USE_PTHREADS
  if(tid > 0)
    {
      int totalLength = 0;

      assert(localTree != tr);

      localTree->manyPartitions          = tr->manyPartitions;
      localTree->NumberOfModels          = tr->NumberOfModels;           
      localTree->rateHetModel            = tr->rateHetModel;
      localTree->saveMemory              = tr->saveMemory;
      localTree->useGappedImplementation = tr->useGappedImplementation;           
      localTree->maxCategories           = tr->maxCategories;      
      localTree->originalCrunchedLength  = tr->originalCrunchedLength;    
      localTree->mxtips                  = tr->mxtips;     
      localTree->numBranches             = tr->numBranches;
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
	    
      localTree->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      localTree->partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));

      if(width > 0 && localTree->saveMemory)
	{
	  localTree->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	    
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
}



int main (int argc, char *argv[])
{ 
  tree  *tr = (tree*)malloc(sizeof(tree));
  
  analdef *adef = (analdef*)malloc(sizeof(analdef));
  
  double **empiricalFrequencies;

  /* not very portable thread to core pinning if PORTABLE_PTHREADS is not defined
     by defualt the cod ebelow is deactivated */

#if (defined(_USE_PTHREADS) && !defined(_PORTABLE_PTHREADS))  
  pinToCore(0);
#endif 

  /* 
     tell the CPU to ignore exceptions generated by denormalized floating point values.
     If this is not done, depending on the input data, the likelihood functions can exhibit 
     substantial run-time differences for vectors of equal length.
  */
     
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
   _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
#endif 

   /* 
      MPI initialization stuff, worry about this later 
   */

#ifdef _FINE_GRAIN_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  printf("\nThis is RAxML FINE-GRAIN MPI Process Number: %d\n", processID);   
  MPI_Barrier(MPI_COMM_WORLD);
#else
  tr->threadID = 0;
#endif


  /* get the start time */

  masterTime = gettime();         
    
  /* initialize the analysis parameters in struct adef to default values */
    
  initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */
  
  get_args(argc, argv, adef, tr); 
  
  /* generate the RAxML output file names and store them in strings */
  
  makeFileNames();
  
  {
    size_t 
      i,
      model;
    
    unsigned char *y;
    
    FILE 
      *byteFile = myfopen(byteFileName, "rb");	 
        
    myBinFread(&(tr->mxtips),                 sizeof(int), 1, byteFile);
    myBinFread(&(tr->originalCrunchedLength), sizeof(int), 1, byteFile);
    myBinFread(&(tr->NumberOfModels),         sizeof(int), 1, byteFile);
    myBinFread(&(tr->gapyness),            sizeof(double), 1, byteFile);
    
    if(adef->perGeneBranchLengths)
      tr->numBranches = tr->NumberOfModels;
    else
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
    
    empiricalFrequencies = (double **)malloc(sizeof(double *) * (size_t)tr->NumberOfModels);
    
    y = (unsigned char *)malloc(sizeof(unsigned char) * ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));
    
    tr->yVector = (unsigned char **)malloc(sizeof(unsigned char *) * ((size_t)(tr->mxtips + 1)));
    
    for(i = 1; i <= (size_t)tr->mxtips; i++)
      tr->yVector[i] = &y[(i - 1) *  (size_t)tr->originalCrunchedLength];	
    
    setupTree(tr);

    /* data structures for convergence criterion need to be initialized after! setupTree */

    if(tr->searchConvergenceCriterion)
      {                     
	tr->bitVectors = initBitVector(tr->mxtips, &(tr->vLength));
	tr->h = initHashTable(tr->mxtips * 4);        
      }

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
    
    myBinFread(y, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips), byteFile);

    fclose(byteFile);
  }
                           


#ifdef _USE_PTHREADS
  /* 
     this main function is the master thread, so if we want to run RAxML with n threads,
     we use startPthreads to start the n-1 worker threads */
  
  startPthreads(tr);
  
  /* via masterBarrier() we invoke parallel regions in which all Pthreads work on computing something, mostly likelihood 
     computations. Have a look at execFunction() in axml.c where we siwtch of the different types of parallel regions.
     
     Although not necessary, below we copy the info stored on tr->partitionData to corresponding copies in each thread.
     While this is shared memory and we don't really need to copy stuff, it was implemented like this to allow for an easier 
     transition to a distributed memory implementation (MPI).
  */
  
  masterBarrier(THREAD_INIT_PARTITION, tr);       
#else      
  /* 
     allocate the required data structures for storing likelihood vectors etc 
  */
  
  initializePartitions(tr, tr, 0, 0);
#endif
       
  /* print out some info on partitions, models, data types etc, not very interesting */
  
  printModelAndProgramInfo(tr, adef, argc, argv);
  
  /* Tells us if the SEV-based memory saving option has been activated in the command line or not.
     printBothOpen() allows to simultaneously print to terminal and to the RAxML_info file, thereby making 
     the terminal output and the info in the RAxML_info file consistent */
  
  printBothOpen("Memory Saving Option: %s\n", (tr->saveMemory == TRUE)?"ENABLED":"DISABLED");   	             
  
  /* 
     initialize model parameters like empirical base frequencies, the rates in the Q matrix, the alpha shape parameters,
     the per-site substitution rates to default starting values */
  
  initModel(tr, empiricalFrequencies);                      
  
  
  
  
  /* 
     this will re-start RAxML exactly where it has left off from a checkpoint file,
     while checkpointing is important and has to be implemented for the library we should not worry about this right now 
  */
  
  if(adef->useCheckpoint)
    {
#ifdef _BAYESIAN
      assert(0);
#endif
      
      /* read checkpoint file */
      restart(tr);
      
      /* continue tree search where we left it off */
      computeBIGRAPID(tr, adef, TRUE); 
    }
  else
    {
      /* not important, only used to keep track of total accumulated exec time 
	 when checkpointing and restarts were used */
      
      accumulatedTime = 0.0;
      
      /* get the starting tree: here we just parse the tree passed via the command line 
	 and do an initial likelihood computation traversal 
	 which we maybe should skeip, TODO */
      
      
      switch(tr->startingTree)
	{
	case randomTree:
	  makeRandomTree(tr);
	  break;
	case givenTree:
	  getStartingTree(tr);     
	  break;
	case parsimonyTree:	     
	  /* runs only on process/thread 0 ! */
	  allocateParsimonyDataStructures(tr);
	  makeParsimonyTreeFast(tr);
	  freeParsimonyDataStructures(tr);
	  break;
	default:
	  assert(0);
	}
      
      /* 
	 here we do an initial full tree traversal on the starting tree using the Felsenstein pruning algorithm 
	 This should basically be the first call to the library that actually computes something :-)
      */
      
      evaluateGeneric(tr, tr->start, TRUE);	 
      
      /* the treeEvaluate() function repeatedly iterates over the entire tree to optimize branch lengths until convergence */
      
      treeEvaluate(tr, 1); 	 	 	 	 	 
      
      /* now start the ML search algorithm */
      
#ifdef _BAYESIAN 
      if(adef->bayesian)
	{
	  /* allocate parsimony data structures for parsimony-biased SPRs */
	  
	  allocateParsimonyDataStructures(tr);
	  mcmc(tr, adef);
	  freeParsimonyDataStructures(tr);
	}
      else
#endif
	computeBIGRAPID(tr, adef, TRUE); 	     
    }            
      
  /* print som more nonsense into the RAxML_info file */
  
  finalizeInfoFile(tr, adef);

  /* return 0 which means that our unix program terminated correctly, the return value is not 1 here */

  return 0;
}


