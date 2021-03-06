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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with 
 *  thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#ifdef __AVX

#ifdef __SIM_SSE3

#define _SSE3_WAS_DEFINED

#undef __SIM_SSE3

#endif

#endif


#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>

#endif


#include "axml.h"



extern const unsigned int mask32[32]; 
/* vector-specific stuff */

extern char **globalArgv;
extern int globalArgc;

#ifdef __SIM_SSE3

#define INTS_PER_VECTOR 4
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128
#define BYTE_ALIGNMENT 16

#endif

#ifdef __AVX

#define INTS_PER_VECTOR 8
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd
#define BYTE_ALIGNMENT 32

#endif

#if !(defined(__SIM_SSE3) || defined(__AVX))
#define BYTE_ALIGNMENT 16
#endif

extern double masterTime;
extern char  workdir[1024];
extern char run_id[128];

/********************************DNA FUNCTIONS *****************************************************************/


static int checkerPars(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checkerPars(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checkerPars(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}

static boolean tipHomogeneityCheckerPars(tree *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->constraintVector[p->number] != grouping) 
	return FALSE;
      else 
	return TRUE;
    }
  else
    {   
      return  (tipHomogeneityCheckerPars(tr, p->next->back, grouping) && tipHomogeneityCheckerPars(tr, p->next->next->back,grouping));      
    }
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->xPars || (s = s->next)->xPars)
    {
      p->xPars = s->xPars;
      s->xPars = 0;
    }

  assert(p->next->xPars || p->next->next->xPars || p->xPars);

}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, boolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->xPars)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
	 computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->xPars) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->xPars) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}



#if (defined(__SIM_SSE3) || defined(__AVX))
#define BIT_COUNT(x, y)  precomputed16_bitcount(x, y)

/* 
   The critical speed of this function
   is mostly a machine/architecure issue
   Nontheless, I decided not to use 
   __builtin_popcount(x)
   for better general portability, albeit it's worth testing
   on x86 64 bit architectures
*/
#else
#define BIT_COUNT(x, y)  precomputed16_bitcount(x, y)
#endif



#if (defined(__SIM_SSE3) || defined(__AVX))

#ifdef __SIM_SSE3

static unsigned int vectorCount(__m128i b)
{
  const unsigned int 
    mu1 = 0x55555555,
    mu2 = 0x33333333,
    mu3 = 0x0F0F0F0F,
    mu4 = 0x0000003F;
  
  unsigned int 
    tcnt[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  
  __m128i 
    m1 = _mm_set_epi32 (mu1, mu1, mu1, mu1),
    m2 = _mm_set_epi32 (mu2, mu2, mu2, mu2),
    m3 = _mm_set_epi32 (mu3, mu3, mu3, mu3),
    m4 = _mm_set_epi32 (mu4, mu4, mu4, mu4),
    tmp1, 
    tmp2;
 

  /* b = (b & 0x55555555) + (b >> 1 & 0x55555555); */
  tmp1 = _mm_srli_epi32(b, 1);                    /* tmp1 = (b >> 1 & 0x55555555)*/
  tmp1 = _mm_and_si128(tmp1, m1); 
  tmp2 = _mm_and_si128(b, m1);                    /* tmp2 = (b & 0x55555555) */
  b    = _mm_add_epi32(tmp1, tmp2);               /*  b = tmp1 + tmp2 */

  /* b = (b & 0x33333333) + (b >> 2 & 0x33333333); */
  tmp1 = _mm_srli_epi32(b, 2);                    /* (b >> 2 & 0x33333333) */
  tmp1 = _mm_and_si128(tmp1, m2); 
  tmp2 = _mm_and_si128(b, m2);                    /* (b & 0x33333333) */
  b    = _mm_add_epi32(tmp1, tmp2);               /* b = tmp1 + tmp2 */

  /* b = (b + (b >> 4)) & 0x0F0F0F0F; */
  tmp1 = _mm_srli_epi32(b, 4);                    /* tmp1 = b >> 4 */
  b = _mm_add_epi32(b, tmp1);                     /* b = b + (b >> 4) */
  b = _mm_and_si128(b, m3);                       /*           & 0x0F0F0F0F */

  /* b = b + (b >> 8); */
  tmp1 = _mm_srli_epi32 (b, 8);                   /* tmp1 = b >> 8 */
  b = _mm_add_epi32(b, tmp1);                     /* b = b + (b >> 8) */
  
  /* b = (b + (b >> 16)) & 0x0000003F; */
  tmp1 = _mm_srli_epi32 (b, 16);                  /* b >> 16 */
  b = _mm_add_epi32(b, tmp1);                     /* b + (b >> 16) */
  b = _mm_and_si128(b, m4);                       /* (b >> 16) & 0x0000003F; */
   
  _mm_store_si128((__m128i *)tcnt, b);

  return tcnt[0] + tcnt[1] + tcnt[2] + tcnt[3];
}

#endif

static inline unsigned int populationCount(INT_TYPE v_N)
{
#ifdef __AVX
  {
    unsigned long int
      res[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
    unsigned int a, b;
    
    _mm256_store_pd((double*)res, v_N);
    
    a = __builtin_popcountl(res[0]) + __builtin_popcountl(res[1]);
    b = __builtin_popcountl(res[2]) + __builtin_popcountl(res[3]);
    
    return (a + b);	   
  }
#else	  
  return (vectorCount(v_N)); 
#endif
}

static void newviewParsimonyIterativeFast(tree *tr)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;

  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];

		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C,
		      v_A, v_C;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);		  		  		  		  
		    
		    v_N = VECTOR_BIT_OR(l_A, l_C);
		    
		    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);		  
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *this[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C, l_G, l_T,
		      v_A, v_C, v_G, v_T;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[2][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[2][i]));
		    l_G = VECTOR_BIT_AND(s_l, s_r);
		    v_G = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[3][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[3][i]));
		    l_T = VECTOR_BIT_AND(s_l, s_r);
		    v_T = VECTOR_BIT_OR(s_l, s_r);
		    
		    v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));	  	 	    	  
		    
		    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
		    VECTOR_STORE((CAST)(&this[2][i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
		    VECTOR_STORE((CAST)(&this[3][i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);	
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[20], 
		      v_A[20];	    	 
		    
		    for(j = 0; j < 20; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < 20; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }
	      }
	      break;
	    default:
	      {
		parsimonyNumber
		  *left[32], 
		  *right[32],
		  *this[32];

		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[32], 
		      v_A[32];	    	 
		    
		    for(j = 0; j < states; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < states; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }	  			
	      }
	    }	  	 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}

static inline unsigned int evaluatePopcount(INT_TYPE v_N, char *precomputed)
{
#ifdef __AVX            	       	   	      
  unsigned long int
    res[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	     
  unsigned int a, b;
	     
  _mm256_store_pd((double*)res, v_N);
  
  a = __builtin_popcountl(res[0]) + __builtin_popcountl(res[1]);
  b = __builtin_popcountl(res[2]) + __builtin_popcountl(res[3]);
	     
  return (a + b);	            
#else      	       
  unsigned int
    sum = 0,
    counts[INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  VECTOR_STORE((CAST)counts, v_N);

  sum += BIT_COUNT(counts[0], precomputed) + BIT_COUNT(counts[1], precomputed);
  sum += BIT_COUNT(counts[2], precomputed) + BIT_COUNT(counts[3], precomputed);          

  return sum;
#endif
}

static unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                       
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),		 
		   v_N = VECTOR_BIT_OR(l_A, l_C);
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N, tr->bits_in_16bits);
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                        
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),
		   l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[2][i])), VECTOR_LOAD((CAST)(&right[2][i]))),
		   l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[3][i])), VECTOR_LOAD((CAST)(&right[3][i]))),
		   v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N, tr->bits_in_16bits);
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i += INTS_PER_VECTOR)
		{                	       
		  int 
		    j;
		  
		  INT_TYPE      
		    l_A,
		    v_N = SET_ALL_BITS_ZERO;     
		  
		  for(j = 0; j < 20; j++)
		    {
		      l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		      v_N = VECTOR_BIT_OR(l_A, v_N);
		    }
		  
		  v_N = VECTOR_AND_NOT(v_N, allOne);
		  
		  sum += evaluatePopcount(v_N, tr->bits_in_16bits);	       
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       *left[32],  
	       *right[32]; 

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	       
		 size_t
		   j;
		 
		 INT_TYPE      
		   l_A,
		   v_N = SET_ALL_BITS_ZERO;     
		 
		 for(j = 0; j < states; j++)
		   {
		     l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		     v_N = VECTOR_BIT_OR(l_A, v_N);
		   }
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N, tr->bits_in_16bits);	       
		 
		 if(sum >= bestScore)	      
		   return sum;		       
	       }
	   }
	 }
    }
  
  return sum;
}


#else

static void newviewParsimonyIterativeFast(tree *tr)
{    
  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];
		
		parsimonyNumber
		   o_A,
		   o_C,
		   t_A,
		   t_C,	
		   t_N;
		
		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];		   

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		  
		    t_N = ~(t_A | t_C);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);		   
		    
		    totalScore += BIT_COUNT(t_N, tr->bits_in_16bits);   
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *this[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		parsimonyNumber
		   o_A,
		   o_C,
		   o_G,
		   o_T,
		   t_A,
		   t_C,
		   t_G,
		   t_T,	
		   t_N;

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];
		    t_G = left[2][i] & right[2][i];	  
		    t_T = left[3][i] & right[3][i];

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		    o_G = left[2][i] | right[2][i];	  
		    o_T = left[3][i] | right[3][i];

		    t_N = ~(t_A | t_C | t_G | t_T);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);
		    this[2][i] = t_G | (t_N & o_G);
		    this[3][i] = t_T | (t_N & o_T); 
		    
		    totalScore += BIT_COUNT(t_N, tr->bits_in_16bits);   
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		parsimonyNumber
		  o_A[20],
		  t_A[20],	  
		  t_N;

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    size_t k;
		    
		    t_N = 0;

		    for(k = 0; k < 20; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;

		    for(k = 0; k < 20; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += BIT_COUNT(t_N, tr->bits_in_16bits); 
		  }
	      }
	      break;
	    default:
	      {		
		parsimonyNumber
		  *left[32],
		  *right[32],
		  *this[32];
		
		parsimonyNumber
		  o_A[32],
		  t_A[32],	  
		  t_N;
		
		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }
		
		for(i = 0; i < width; i++)
		  {	 	  
		    t_N = 0;
		    
		    for(k = 0; k < states; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;
		    
		    for(k = 0; k < states; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += BIT_COUNT(t_N, tr->bits_in_16bits); 
		  }
	      }			      
	    } 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



static unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber 
	       t_A,
	       t_C,	      
	       t_N,
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i++)
	       {                	                       
		 t_A = left[0][i] & right[0][i];
		 t_C = left[1][i] & right[1][i];
		 
		  t_N = ~(t_A | t_C);

		  sum += BIT_COUNT(t_N, tr->bits_in_16bits);    
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       t_A,
	       t_C,
	       t_G,
	       t_T,
	       t_N,
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i++)
	       {                	                        
		  t_A = left[0][i] & right[0][i];
		  t_C = left[1][i] & right[1][i];
		  t_G = left[2][i] & right[2][i];	  
		  t_T = left[3][i] & right[3][i];

		  t_N = ~(t_A | t_C | t_G | t_T);

		  sum += BIT_COUNT(t_N, tr->bits_in_16bits);     
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i++)
		{ 
		  t_N = 0;
		  
		  for(k = 0; k < 20; k++)
		    {
		      t_A = left[k][i] & right[k][i];
		      t_N = t_N | t_A;
		    }
  	       
		  t_N = ~t_N;

		  sum += BIT_COUNT(t_N, tr->bits_in_16bits);      
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[32], 
	       *right[32];  

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i++)
	       {                	       
		 t_N = 0;
		  
		 for(k = 0; k < states; k++)
		   {
		     t_A = left[k][i] & right[k][i];
		     t_N = t_N | t_A;
		   }
  	       
		  t_N = ~t_N;

		  sum += BIT_COUNT(t_N, tr->bits_in_16bits);      
		  		  		 
		 if(sum >= bestScore)			  
		   return sum;			   
	       }	     	     
	   }
	 }
    }
  
  return sum;
}

#endif






static unsigned int evaluateParsimony(tree *tr, nodeptr p, boolean full)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->xPars)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->xPars)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;

  result = evaluateParsimonyIterativeFast(tr);

  return result;
}


static void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr);      
  }
}





/****************************************************************************************************************************************/

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches); 
   
  newviewParsimony(tr, p);     
} 

/*
  static nodeptr buildNewTip (tree *tr, nodeptr p)
  { 
  nodeptr  q;
  
  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
  assert(q == q->next->next->next);
  assert(q->xPars || q->next->xPars || q->next->next->xPars);
  return  q;
  } 
*/

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}


static void testInsertParsimony (tree *tr, nodeptr p, nodeptr q, boolean saveBranches)
{ 
  unsigned int 
    mp;
 
  nodeptr  
    r = q->back;   

  boolean 
    doIt = TRUE;
    
  if(tr->grouped)
    {
      int 
	rNumber = tr->constraintVector[r->number],
	qNumber = tr->constraintVector[q->number],
	pNumber = tr->constraintVector[p->number];

      doIt = FALSE;
     
      if(pNumber == -9)
	pNumber = checkerPars(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checkerPars(tr, q);

	  if(rNumber == -9)
	    rNumber = checkerPars(tr, r);

	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;       
	}
    }

  if(doIt)
    {
      double 
	z[NUM_BRANCHES];
      
      if(saveBranches)
	{
	  int i;
	  
	  for(i = 0; i < tr->numBranches; i++)
	    z[i] = q->z[i];
	}

      insertParsimony(tr, p, q);   
  
      mp = evaluateParsimony(tr, p->next->next, FALSE);                      

      if(mp < tr->bestParsimony)
	{
	  tr->bestParsimony = mp;
	  tr->insertNode = q;
	  tr->removeNode = p;
	}
      
      if(saveBranches)
	hookup(q, r, z, tr->numBranches);
      else
	hookupDefault(q, r, tr->numBranches);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
    }
       
  return;
} 


static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr); 
}


static void addTraverseParsimony (tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav, boolean doAll, boolean saveBranches)
{        
  if (doAll || (--mintrav <= 0))               
    testInsertParsimony(tr, p, q, saveBranches);	                 

  if (((q->number > tr->mxtips)) && ((--maxtrav > 0) || doAll))
    {	      
      addTraverseParsimony(tr, p, q->next->back, mintrav, maxtrav, doAll, saveBranches);	      
      addTraverseParsimony(tr, p, q->next->next->back, mintrav, maxtrav, doAll, saveBranches);              	     
    }
}




static void makePermutationFast(int *perm, int n, tree *tr)
{    
  int  
    i, 
    j, 
    k;

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum(&tr->randomNumberSeed);

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

static nodeptr  removeNodeParsimony (nodeptr p, tree *tr)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupDefault(q, r, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}

static int rearrangeParsimony(tree *tr, nodeptr p, int mintrav, int maxtrav, boolean doAll)  
{   
  nodeptr  
    p1, 
    p2, 
    q, 
    q1, 
    q2;
  
  int      
    mintrav2; 

  boolean 
    doP = TRUE,
    doQ = TRUE;
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  assert(mintrav == 1);

  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(tr->constrained)
    {    
      if(! tipHomogeneityCheckerPars(tr, p->back, 0))
	doP = FALSE;
	
      if(! tipHomogeneityCheckerPars(tr, q->back, 0))
	doQ = FALSE;
		        
      if(doQ == FALSE && doP == FALSE)
	return 0;
    }  

  if((p->number > tr->mxtips) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if ((p1->number > tr->mxtips) || (p2->number > tr->mxtips)) 
	{	  	  
	  removeNodeParsimony(p, tr);	  	 

	  if ((p1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p1->next->back, mintrav, maxtrav, doAll, FALSE);         
	      addTraverseParsimony(tr, p, p1->next->next->back, mintrav, maxtrav, doAll, FALSE);          
	    }
	 
	  if ((p2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p2->next->back, mintrav, maxtrav, doAll, FALSE);
	      addTraverseParsimony(tr, p, p2->next->next->back, mintrav, maxtrav, doAll, FALSE);          
	    }
	    
	   
	  hookupDefault(p->next,       p1, tr->numBranches); 
	  hookupDefault(p->next->next, p2, tr->numBranches);	   	    	    

	  newviewParsimony(tr, p);
	}
    }  
       
  if ((q->number > tr->mxtips) && (maxtrav > 0) && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

      if (
	  (
	   (q1->number > tr->mxtips) && 
	   ((q1->next->back->number > tr->mxtips) || (q1->next->next->back->number > tr->mxtips))
	   )
	  ||
	  (
	   (q2->number > tr->mxtips) && 
	   ((q2->next->back->number > tr->mxtips) || (q2->next->next->back->number > tr->mxtips))
	   )
	  )
	{	   

	  removeNodeParsimony(q, tr);
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if ((q1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q1->next->back, mintrav2 , maxtrav, doAll, FALSE);
	      addTraverseParsimony(tr, q, q1->next->next->back, mintrav2 , maxtrav, doAll, FALSE);         
	    }
	 
	  if ((q2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q2->next->back, mintrav2 , maxtrav, doAll, FALSE);
	      addTraverseParsimony(tr, q, q2->next->next->back, mintrav2 , maxtrav, doAll, FALSE);          
	    }	   
	   
	  hookupDefault(q->next,       q1, tr->numBranches); 
	  hookupDefault(q->next->next, q2, tr->numBranches);
	   
	  newviewParsimony(tr, q);
	}
    }

  return 1;
} 


static void restoreTreeRearrangeParsimony(tree *tr)
{    
  removeNodeParsimony(tr->removeNode, tr);  
  restoreTreeParsimony(tr, tr->removeNode, tr->insertNode);  
}

/*
static boolean isInformative2(tree *tr, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = 15;

  unsigned char
    nucleotide,
    target = 0;
  	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;      	           
    }
  
  
  if(check[1] > 1)
    {
      informativeCounter++;    
      target = target | 1;
    }
  if(check[2] > 1)
    {
      informativeCounter++; 
      target = target | 2;
    }
  if(check[4] > 1)
    {
      informativeCounter++; 
      target = target | 4;
    }
  if(check[8] > 1)
    {
      informativeCounter++; 
      target = target | 8;
    }
	  
  if(informativeCounter >= 2)
    return TRUE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(j == 3 || j == 5 || j == 6 || j == 7 || j == 9 || j == 10 || j == 11 || 
	     j == 12 || j == 13 || j == 14)
	    {
	      if(check[j] > 1)
		{
		  if(!(target & j))
		    return TRUE;
		}
	    }
	} 
    }
     
  return FALSE;	     
}
*/

static boolean isInformative(tree *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}


static void determineUninformativeSites(tree *tr, int *informative)
{
  int 
    model,
    number = 0,
    i;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */


  for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	   if(isInformative(tr, tr->partitionData[model].dataType, i))
	     informative[i] = 1;
	   else
	     {
	       informative[i] = 0;
	       number++;
	     }  
	}      
    }

 
 
  /* printf("Uninformative Patterns: %d\n", number); */
}


static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if((p->number <= tr->mxtips))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	}            
     
      assert(found != 0);

      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}



static void nodeRectifierPars(tree *tr)
{
  nodeptr *np = (nodeptr *)malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = FALSE;

  /* TODO why is tr->rooted set to FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 

 
  rax_free(np);
}


  
static void compressDNA(tree *tr, int *informative)
{
  size_t
    totalNodes,
    i,
    model;
   
  totalNodes = 2 * (size_t)tr->mxtips;

 

  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = (size_t)tr->partitionData[model].states,       
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)malloc(states * sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)malloc(states * sizeof(parsimonyNumber));
      
      for(i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->aliaswgt[i];     
  
      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SIM_SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      
      tr->partitionData[model].parsVect = (parsimonyNumber *)malloc_aligned((size_t)compressedEntriesPadded * states * totalNodes * sizeof(parsimonyNumber));
     
      for(i = 0; i < compressedEntriesPadded * states * totalNodes; i++)      
	tr->partitionData[model].parsVect[i] = 0;          

      for(i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(tr->partitionData[model].parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const unsigned int 
		    *bitValue = getBitVector(tr->partitionData[model].dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & mask32[k])
			    compressedValues[k] |= mask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= mask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	}               
  
      tr->partitionData[model].parsimonyLength = compressedEntriesPadded;   

      rax_free(compressedTips);
      rax_free(compressedValues);
    }
  
  tr->parsimonyScore = (unsigned int*)malloc_aligned(sizeof(unsigned int) * totalNodes);  
          
  for(i = 0; i < totalNodes; i++) 
    tr->parsimonyScore[i] = 0;
}



static void stepwiseAddition(tree *tr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;
    
  mp = evaluateParsimonyIterativeFast(tr);
  
  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(tr, p, q->next->back);	      
      stepwiseAddition(tr, p, q->next->next->back);              	     
    }
}




void allocateParsimonyDataStructures(tree *tr)
{
  int 
    i,
    *informative = (int *)malloc(sizeof(int) * (size_t)tr->originalCrunchedLength);
 
  determineUninformativeSites(tr, informative);     

  compressDNA(tr, informative);

  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 1; i++)
    {
      nodeptr 
	p = tr->nodep[i];

      p->xPars = 1;
      p->next->xPars = 0;
      p->next->next->xPars = 0;
    }

  tr->ti = (int*)malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  rax_free(informative); 
}

void pllFreeParsimonyDataStructures(tree *tr)
{
  size_t 
    model;

  rax_free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    rax_free(tr->partitionData[model].parsVect);
  
  rax_free(tr->ti);
}


void pllMakeParsimonyTreeFast(tree *tr)
{   
  nodeptr  
    p, 
    f;    

  int 
    i, 
    nextsp,
    *perm        = (int *)malloc((size_t)(tr->mxtips + 1) * sizeof(int));  

  unsigned int 
    randomMP, 
    startMP;         
  
  assert(!tr->constrained);

  makePermutationFast(perm, tr->mxtips, tr);
  
  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
  
  buildSimpleTree(tr, perm[1], perm[2], perm[3]);      
  
  f = tr->start;       
  
  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;
      
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped)
	{
	  int 
	    number = p->back->number;	  	 

	  tr->constraintVector[number] = -9;
	}
          
      stepwiseAddition(tr, q, f->back);      	  	 
      
      {
	nodeptr	  
	  r = tr->insertNode->back;
	
	int counter = 4;
	
	hookupDefault(q->next,       tr->insertNode, tr->numBranches);
	hookupDefault(q->next->next, r, tr->numBranches);
	
	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(tr);	
      }
    }    
  
  printf("ADD: %d\n", tr->bestParsimony); 
  
  nodeRectifierPars(tr);
  
  randomMP = tr->bestParsimony;        
  
  do
    {
      startMP = randomMP;
      nodeRectifierPars(tr);
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{
	  rearrangeParsimony(tr, tr->nodep[i], 1, 20, FALSE);
	  if(tr->bestParsimony < randomMP)
	    {		
	      restoreTreeRearrangeParsimony(tr);
	      randomMP = tr->bestParsimony;
	    }
	}      		  	   
    }
  while(randomMP < startMP);
  
  printf("OPT: %d\n", tr->bestParsimony);
} 

void parsimonySPR(nodeptr p, tree *tr)
{
  int i;

  double   
    p1z[NUM_BRANCHES], 
    p2z[NUM_BRANCHES];

  nodeptr 
    p1 = p->next->back,
    p2 = p->next->next->back;

  unsigned int score = evaluateParsimony(tr, p, TRUE);

  printf("parsimonyScore: %u\n", score);

  for(i = 0; i < tr->numBranches; i++)
    {
      p1z[i] = p1->z[i];
      p2z[i] = p2->z[i];	   	   
    }
  
  tr->bestParsimony = INT_MAX; 

  hookupDefault(p1, p2, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;

  if (p1->number > tr->mxtips) 
    {
      addTraverseParsimony(tr, p, p1->next->back, 0, 0, TRUE, TRUE);         
      addTraverseParsimony(tr, p, p1->next->next->back, 0, 0, TRUE, TRUE);          
    }
  
  if(p2->number > tr->mxtips)
    {
      addTraverseParsimony(tr, p, p2->next->back, 0, 0, TRUE, TRUE);
      addTraverseParsimony(tr, p, p2->next->next->back, 0, 0, TRUE, TRUE);          
    }

  printf("best %u nodes %d %d\n",tr->bestParsimony, tr->insertNode->number, tr->insertNode->back->number);

  hookup(p1, p->next, p1z,       tr->numBranches);
  hookup(p2, p->next->next, p2z, tr->numBranches);
}


#ifdef __AVX

#ifdef _SSE3_WAS_DEFINED

#define __SIM_SSE3

#undef _SSE3_WAS_DEFINED

#endif

#endif
