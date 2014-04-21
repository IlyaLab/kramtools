
#ifndef _analysis_h_
#define _analysis_h_

#ifdef __cplusplus
extern "C" {
#endif

struct CovariateAnalysis {

#define FAIL_SAMPLES 0x00000001 // Too few samples after filtering.
#define FAIL_E_DEGEN 0x00000002 // Early degeneracy, before filtering.
#define FAIL_L_DEGEN 0x00000004 // Late degeneracy, after filtering.
#define FAIL_TOOMANY 0x00000008 // Too many categories
#define FAIL_MATH    0x00000010 // Something went south in the math.

	unsigned status;

	/**
	  * These are defined as bitfields by mtm library header (mttypeid.h).
	  * Yes, there's a reason for bitfields, not enum!
	  * See documentation from mtm library.
	  */
	struct {
		unsigned left, right;
	} stat_class;

	/**
	  * This array deals with the "wasted" samples from the two covariates.
	  * If '+' means present, '-' missing...
	  *
	  * covariate 1: - - + - - - + + + - + + - - + + - -
	  * covariate 2: - + + + + - + + + - + + - - - - - - 
	  *                  1       1 1 1   1 1     0 0      unused == 2
	  *                0 1 0 0   1 1 1   1 1              unused == 3
	  */

	struct {
		signed unused; // signed to contain invalid vals to guarantee init.
		struct Statistic result;
	} waste[2];

	/**
	  * Indicates the sign of the relationship if both variables are 
	  * ordinal: Con-Con, Con-Ord, or Con-Cat where Cat is boolean
	  * (so an order is implied). Unused if either feature is categorical
	  * with >2 categories.
	  */
	int sign;

	struct Statistic result;
};

/**
  * This pre-allocates all the working memory that will be needed.
  */
int  covan_init( int columns );

/**
 * Returns non-zero on error and 0 otherwise.
 */
int  covan_exec( const struct mt_row_pair *pair, struct CovariateAnalysis * );

#ifdef __cplusplus
}
#endif

#endif

