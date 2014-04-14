
#ifndef _analysis_h_
#define _analysis_h_

#ifdef __cplusplus
extern "C" {
#endif


/**
 * This is intended to provide a statistic-independent way to return
 * (dis)qualifying information about the samples for which a statistic
 * was calculated.
 */

/**
  * Left triple identifies feature 1 and right triple feature 2, so"CatNum"
  * means feature 1 (also left in method call argument order) is categorical 
  * and feature 2 (right in method calls) is numeric.
  */
enum CovarTypes {
	UnknownCovar = 0,
	NumNum,     // numeric,     numeric
	CatNum,     // categorical, numeric
	NumCat,     // numeric,     categorical 
	CatCat,     // categorical, categorical
	IndetDegen, // indeterminate, degenerate
	DegenIndet, // degenerate,  indeterminate 
	CovarTypeCount
};

extern const char *COVAR_TYPE_STR[];

struct CommonStats {

	/**
	  * Number of covariate pairs in which BOTH values were non-missing.
	  * ...signed so that it can contain a clearly invalid value to
	  * insure proper initialization.
	  */
	signed N;

	/**
	  * P-value of whatever test was performed.
	  */
	double P;
};
#define CLEAR_COMMON(p) {(p)->N=-1;(p)->P=nan("nan");}

enum HypothesisTest {
	UnknownHypTest = 0,
	FisherExact,
	ChiSquare,
	MannWhitney,
	KruskalWallis,
	Spearman,
	Pearson,
	HypothesisTestCount
};

extern const char *HypTestNames[];

struct SpearmanStats {
	double   rho;
	unsigned ties[2];
};

struct PearsonStats {
	double   rho;
};

struct KruskalWallisStats {
	double   K;
	unsigned ties; // for numeric covariate
};

struct MannWhitneyStats {
	double   U;
	unsigned ties; // for numeric covariate
};

struct ChiSquareStats {
	double   chisq;
	double   min_expected_cell_count;
	unsigned rows; // ...AFTER any merging and/or culling was
	unsigned cols; // performed, i.e. on which P was calculated.
	unsigned empty_cells;
};

struct FisherExactStats {
	unsigned rows;
	unsigned cols;
};


struct Analysis {

	unsigned status;

	struct CommonStats common;

	/**
	  * Spearman correlation -- valid when neither covariate is a 
	  * categorical variable with >2 categories. That is, if
	  * 1) both covariates are numeric
	  * 2) one covariate is numeric and one is a binary categorical
	  * 3) both covariates are binary categorical.
	  */
	double spearman_rho;

	enum CovarTypes kind;

	/**
	  * What was actually performed to generate the P-value in common 
	  * above. (There is no P-value associated with spearman_rho except
	  * when it coincides with spearman.rho below.)
	  */
	enum HypothesisTest test;

	union {
		struct SpearmanStats      spearman;
		struct PearsonStats       pearson;
		struct KruskalWallisStats kruskal;
		struct ChiSquareStats     chisq;
		struct FisherExactStats   fisherx;
	};

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
		struct CommonStats common;
		struct KruskalWallisStats kruskal;
	} waste[2];

	/**
	 * catcovars.cpp uses this to return information about row/column 
	 * culling and/or merging. No other uses (currently).
	 */
#define SUMMARY_LOG_LEN  (256)

	char log[ SUMMARY_LOG_LEN ];
};

void clear_summary( struct Analysis *cs );

#define EOLOG         ('-')

int analysis_init( int columns );

#define FAIL_SAMPLES 0x00000001 // Too few samples after filtering.
#define FAIL_E_DEGEN 0x00000002 // Early degeneracy, before filtering.
#define FAIL_L_DEGEN 0x00000004 // Late degeneracy, after filtering.
#define FAIL_TOOMANY 0x00000008 // Too many categories
#define FAIL_MATH    0x00000010 // Something went south in the math.

/**
 * Returns non-zero on error and 0 otherwise.
 */
unsigned analyze_pair( const struct mt_row_pair *pair,
		struct Analysis * );

#ifdef __cplusplus
}
#endif

#endif

