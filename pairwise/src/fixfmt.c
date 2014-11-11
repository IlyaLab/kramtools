
/**
  * This module is deprecated and should be removed eventually.
  * It's capabilities are subsumed by the varfmt.c module.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "mtmatrix.h"
#include "featpair.h"
#include "mtsclass.h"
#include "stattest.h"
#include "analysis.h"
#include "fixfmt.h"
#include "varfmt.h" // even these legacy methods must follow the new EMITTER_SIG
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
static double __attribute__((always_inline)) _clamped_neglog( double prob ) {
	prob = -log10( prob );
	return isinf(prob) || (prob > 300 ) ? 300.0 : prob;
}


/**
  * This need not support all possibilities since this is transitional
  * code.
  * analysis_type(NN,CN,NC,CC,1X,X2):hypothesis_test
  */
static const char *COVAR_TYPE_STR( unsigned l, unsigned r ) {
	if( MTM_STATCLASS_CATEGORICAL == l ) {
		if( MTM_STATCLASS_CATEGORICAL == r ) {
			return "CC";
		} else {
			assert( MTM_STATCLASS_CONTINUOUS == r );
			return "CN";
		}
	} else {
		assert( MTM_STATCLASS_CONTINUOUS == l );
		if( MTM_STATCLASS_CATEGORICAL == r ) {
			return "NC";
		} else {
			assert( MTM_STATCLASS_CONTINUOUS == r );
			return "NN";
		}
	}
}


void format_tcga( EMITTER_SIG ) {

	const int u0
		= covan->waste[0].unused;
	const int u1
		= covan->waste[1].unused;

	if( pair->l.name != NULL && pair->r.name != NULL ) {
		fprintf( fp, "%s\t%s\t",  pair->l.name,   pair->r.name );
	} else {
		fprintf( fp, "%d\t%d\t",  pair->l.offset, pair->r.offset );
	}

	fprintf( fp, 
		"%s\t"
		"%+.2f\t"
		"%d\t"
		"%.3f\t"
		"%d\t%.3f\t"
		"%d\t%.3f\t"
		"%s\n", 
		COVAR_TYPE_STR( covan->stat_class.left, covan->stat_class.right ),
		covan->sign, 
		covan->result.sample_count, 
		_clamped_neglog( covan->result.probability ),
		u0, u0 > 0 ? _clamped_neglog( covan->waste[0].result.probability ) : 0.0,
		u1, u1 > 0 ? _clamped_neglog( covan->waste[1].result.probability ) : 0.0,
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
void format_standard( EMITTER_SIG ) {

	if( pair->l.name != NULL && pair->r.name != NULL ) {
		fprintf( fp, "%s\t%s\t",  pair->l.name,   pair->r.name );
	} else {
		fprintf( fp, "%d\t%d\t",  pair->l.offset, pair->r.offset);
	}

	fprintf( fp, 
		"%s:%s\t"
		"%+.2f\t"
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


void format_abbreviated( EMITTER_SIG ) {

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

