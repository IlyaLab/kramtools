
#ifndef _analysis_h_
#define _analysis_h_

#ifdef __cplusplus
extern "C" {
#endif


enum FeatureClass {
	Unknown = 0,
	Categorical,
	Ordinal,
	Continuous,
	FeatureClassCount
};

struct CovariateAnalysis {

#define FAIL_SAMPLES 0x00000001 // Too few samples after filtering.
#define FAIL_E_DEGEN 0x00000002 // Early degeneracy, before filtering.
#define FAIL_L_DEGEN 0x00000004 // Late degeneracy, after filtering.
#define FAIL_TOOMANY 0x00000008 // Too many categories
#define FAIL_MATH    0x00000010 // Something went south in the math.

	unsigned status;

	enum FeatureClass lclass, rclass;

	/**
	  * Spearman correlation -- valid when neither covariate is a 
	  * categorical variable with >2 categories. That is, if
	  * 1) both covariates are numeric
	  * 2) one covariate is numeric and one is a binary categorical
	  * 3) both covariates are binary categorical.
	  */
	struct Statistic correlation;

	struct Statistic association;


	/**
	  * This array deals with the "wasted" samples from the two covariates.
	  * If '+' means present, '-' missing...
	  *
	  * covariate 1: - - + - - - + + + - + + - - + + - -
	  * covariate 2: - + + + + - + + + - + + - - - - - - 
	  *                  1       1 1 1   1 1     0 0      unused == 2
	  *                0 1 0 0   1 1 1   1 1              unused == 3
	  */

	struct Statistic waste[2];
/*
	struct {
		signed unused; // signed to contain invalid vals to guarantee init.
		struct CommonStats common;
		struct KruskalWallisStats kruskal;
	} waste[2];
*/
};

void covan_clear( struct CovariateAnalysis *cs );

#define EOLOG         ('-')

int covan_init( int columns );

/**
 * Returns non-zero on error and 0 otherwise.
 */
int covan_exec( const struct mt_row_pair *pair, struct CovariateAnalysis * );

#ifdef __cplusplus
}
#endif

#endif

