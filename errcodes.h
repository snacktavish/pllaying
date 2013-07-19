#ifndef ERRCODES_H
#define ERRCODES_H

#define  PLL_NNI_P_TIP                  1 << 0          /**< Node p is a tip */
#define  PLL_NNI_Q_TIP                  1 << 1          /**< Node p->back is a tip */

#define  PLL_PARTITION_OUT_OF_BOUNDS    1 << 0      /**< Trying to access a partition index that is out of bounds */
#define  PLL_BASE_FREQUENCIES_DO_NOT_SUM_TO_1 1 << 1      /**< base frequencies don't sum to 1.0 */
#endif
