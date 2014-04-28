
#ifndef _analysis_h_
#define _analysis_h_

#ifdef __cplusplus
extern "C" {
#endif

struct CovariateAnalysis {

// In approximate order of likelihood...
#define COVAN_E_SAMPLES_SIZE 0x00000001 // Too few samples
#define COVAN_E_UNIVAR_DEGEN 0x00000002 // Univariate degeneracy
#define COVAN_E_COVAR_DEGEN  0x00000004 // Covariate degeneracy
#define COVAN_E_MATH         0x00000008 // Something went south in the math.
#define COVAN_E_TOOMANY_CATS 0x00000010 // Too many categories
#define COVAN_E_MASK         0x0000001F
// This list must stay synchronized with the text in usage_full.txt.

	unsigned status;

	/**
	  * These are defined as bitfields by mtm library header (mtsclass.h).
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
	float sign;

	struct Statistic result;
};
typedef struct CovariateAnalysis CovariateAnalysis_t;
typedef const CovariateAnalysis_t COVARIATEANALYSIS_T;
/**
  * This pre-allocates all the working memory that will be needed.
  */
int  covan_init( int columns );

/**
 * Returns non-zero on error and 0 otherwise.
 */
int  covan_exec( const struct feature_pair *pair, struct CovariateAnalysis * );

#ifdef __cplusplus
}
#endif

#endif

