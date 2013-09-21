#ifndef __pll_ALIGNMENT__
#define __pll_ALIGNMENT__

/** @brief Generic structure for storing a multiple sequence alignment */
typedef struct
 {
   int              sequenceCount;      /**< @brief Number of sequences */
   int              sequenceLength;     /**< @brief Length of sequences */
   char          ** sequenceLabels;     /**< @brief An array of where the \a i-th element is the name of the \a i-th sequence */
   unsigned char ** sequenceData;       /**< @brief The actual sequence data */
   int            * siteWeights;        /**< @brief An array where the \a i-th element indicates how many times site \a i appeared (prior to duplicates removal) in the alignment */
 } pllAlignmentData;

#endif
