
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "mtmatrix.h"
#include "mttypeid.h"
#include "stattest.h"
#include "analysis.h"
 
/**
  * "Sheila's format"
  * 1  -- feature A
  * 2  -- feature B
  * 3  -- analysis type (NN,CN,NC,CC,1X,X2)
  *      N: numeric
  *      C: categorical
  *      1X (in categorical pairs) feature B only had 1 (surviving) category
  *      X2 (in categorical pairs) feature A only had 1 (surviving) category
  * 4  -- Spearman correlation coefficient
  * 5  -- number of samples used in pairwise test
  * 6  -- -log10(p-value) (of whatever test was performed)
  * 7  -- count of unused samples in feature A
  * 8  -- -log10(p-value) of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature A
  * 9  -- count of unused samples in feature B
  * 10 -- -log10(p-value) of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature B
  */
static double __attribute__((always_inline)) _clamped( double prob ) {
	prob = -log10( prob );
	return isinf(prob) || (prob > 300 ) ? 300.0 : prob;
}


/**
  * This need not support all possibilities since this is transitional
  * code.
  * analysis_type(NN,CN,NC,CC,1X,X2):hypothesis_test
  */
static const char *COVAR_TYPE_STR( unsigned l, unsigned r ) {
	if( STAT_CLASS_CATEGORICAL == l ) {
		if( STAT_CLASS_CATEGORICAL == r ) {
			return "CC";
		} else {
			assert( STAT_CLASS_CONTINUOUS == r );
			return "CN";
		}
	} else {
		assert( STAT_CLASS_CONTINUOUS == l );
		if( STAT_CLASS_CATEGORICAL == r ) {
			return "NC";
		} else {
			assert( STAT_CLASS_CONTINUOUS == r );
			return "NN";
		}
	}
}


void format_tcga( 
		const struct mt_row_pair *pair, 
		const struct CovariateAnalysis *covan, FILE *fp ) {

	double NLOGP[3];
	char rho[8];

	if( covan->sign == 0 )
		strcpy( rho, "NA" );
	else 
		sprintf( rho, "%+d", covan->sign );

	NLOGP[0] = _clamped( covan->result.probability );
	NLOGP[1] = _clamped( covan->waste[0].result.probability );
	NLOGP[2] = _clamped( covan->waste[1].result.probability );

	if( pair->l.name != NULL && pair->r.name != NULL ) {
		fprintf( fp, "%s\t%s\t",  pair->l.name,   pair->r.name );
	} else {
		fprintf( fp, "%d\t%d\t",  pair->l.offset, pair->r.offset );
	}

	fprintf( fp, 
		"%s\t"
		"%s\t"
		"%d\t"
		"%.3f\t"
		"%d\t%.3f\t"
		"%d\t%.3f\t"
		"%s\n", 
		COVAR_TYPE_STR( covan->stat_class.left, covan->stat_class.right ),
		rho, 
		covan->result.sample_count, 
		NLOGP[0       ],
		covan->waste[0].unused, NLOGP[1],
		covan->waste[1].unused, NLOGP[2],
		covan->result.log );
}



/**
  * Standard format
  * 1  -- feature A
  * 2  -- feature B
  * 3  -- analysis_type(NN,CN,NC,CC,1X,X2):hypothesis_test
  *      N: numeric
  *      C: categorical
  *      1X (in categorical pairs) feature B only had 1 (surviving) category
  *      X2 (in categorical pairs) feature A only had 1 (surviving) category
  *      Hypothesis test is a short unambiguous string hyphenated onto the
  *      analysis type so that it's still one column (for now).
  * 4  -- Spearman correlation coefficient
  * 5  -- number of samples used in pairwise test
  * 6  -- p-value
  * 7  -- count of unused samples in feature A
  * 8  -- p-value of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature A
  * 9  -- count of unused samples in feature B
  * 10 -- p-value of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature B
  * 11 -- log (of contingency table culling)
  */
void format_standard( 
		const struct mt_row_pair *pair, 
		const struct CovariateAnalysis *covan, FILE *fp ) {

	if( pair->l.name != NULL && pair->r.name != NULL ) {
		fprintf( fp, "%s\t%s\t",  pair->l.name,   pair->r.name );
	} else {
		fprintf( fp, "%d\t%d\t",  pair->l.offset, pair->r.offset);
	}

	fprintf( fp, 
		"%s:%s\t"
		"%+d\t"
		"%d\t"
		"%.3e\t"
		"%d\t%.3e\t"
		"%d\t%.3e\t"
		"%s\n", 
		COVAR_TYPE_STR( covan->stat_class.left, covan->stat_class.right ), covan->result.name,
		covan->sign, 
		covan->result.sample_count, 
		covan->result.probability,
		covan->waste[0].unused, covan->waste[0].result.probability,
		covan->waste[1].unused, covan->waste[1].result.probability,
		covan->result.log );
}


void format_abbreviated( 
		const struct mt_row_pair *pair, 
		const struct CovariateAnalysis *covan, FILE *fp ) {

	if( pair->l.name != NULL && pair->r.name != NULL )
		fprintf( fp, "%s\t%s\t", pair->l.name, pair->r.name );
	else
		fprintf( fp, "%d\t%d\t",  pair->l.offset, pair->r.offset);

	fprintf( fp, 
		"%s\t"
		"%04x\t"
		"%d\t%.3e\t"
		"%d\t%.3e\t%d\t"
		"%d\t%.3e\t%d\t"
		"%.3f\t"
		"%s\n", 
		COVAR_TYPE_STR( covan->stat_class.left, covan->stat_class.right ),
		covan->status,
		covan->result.sample_count,          covan->result.probability,
		covan->waste[0].result.sample_count, covan->waste[0].result.probability, covan->waste[0].unused,
		covan->waste[1].result.sample_count, covan->waste[1].result.probability, covan->waste[1].unused,
		covan->sign == 0 ? +nan("nan") : covan->sign, 
		covan->result.log );
}

