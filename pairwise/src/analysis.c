/**
 * This module is really just a dispatcher and metadata collector for
 * the four currently supported feature class combinations:
 * 1) continuous, continuous
 * 2) categorical, continuous
 * 3) continuous, categorical 
 * 4) categorical, categorical
 *
 * Importantly (for performance) each (statically-allocated) class 
 * pre-allocates all the buffer space it will require. And nothing
 * is freed until process exit.
 * 
 * Re: performance, it may appear wasteful that I'm actually copying
 * values from the input arrays *into* the class instances; it's
 * certainly not ideal. But a significant part of this code's work
 * is -filtering- out just the covariate pairs that are BOTH present.
 * In this context, copying isn't so silly, especially as recipient
 * class instances structure the copied data on the fly for optimal
 * downstream computation.
 *
 * Code here is motivated by two overriding concerns:
 * 1. the statistical processing classes mix/cat/num should depend
 *    on nothing more than the struct Statistic class.
 * 2. everything on data should be accomplished in 1-pass! This
 *    puts some hard constraints on the approach to code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "stattest.h"
#include "mtmatrix.h"
#include "featpair.h"
#include "analysis.h"

#include "cat.h"
#include "mix.h"
#include "num.h"
#include "args.h"
#include "mtsclass.h"
#include "limits.h"


/**
  * Most arrays are pre-allocated and sized according to this variable.
  */
static int max_sample_count = 0;

/**
 * These classes handle the actual feature1 vs feature2 analyses.
 */
static void *_caccum = NULL;
static void *_maccum = NULL;
static void *_naccum = NULL;

/**
 * These classes handle comparisons between the two groups WITHIN a
 * single feature distinguished according to whether the corresponding
 * covariate (in the other feature) is present (not "NA"). In continuous/
 * continuous comparisons are both used, in continuous/categorical
 * comparisons one or the other is used, and in cat/cat neither is used.
 *
 * Notice that in general counts by category in categorical/continuous
 * pairs are generally unavailable with one exception: the category 0's
 * count is always output by the Kruskal-Wallis test. This is motivated
 * primarily by the application defined above. As a result, it is 
 * ESSENTIAL THAT SAMPLES PRESENT IN FEAT1 BUT MISSING IN FEAT2 ARE
 * TAGGED WITH 0 IN _Lwaste (vica versa for 2).
 */
static void *_Lwaste = NULL;
static void *_Rwaste = NULL;

////////////////////////////////////////////////////////////////////////////
// Public API
////////////////////////////////////////////////////////////////////////////

void covan_fini( void ) {

	if( _Rwaste ) { mix_destroy( _Rwaste ); _Rwaste = NULL; }
	if( _Lwaste ) { mix_destroy( _Lwaste ); _Lwaste = NULL; }
	if( _naccum ) { con_destroy( _naccum ); _naccum = NULL; }
	if( _maccum ) { mix_destroy( _maccum ); _maccum = NULL; }
	if( _caccum ) { cat_destroy( _caccum ); _caccum = NULL; }
}

/**
 * This must pre-allocate all the memory we might need for every combination
 * of data types.
 */
int covan_init( int columns ) {

	max_sample_count = columns; // EVERYTHING depends on this.

	_caccum = cat_create( MAX_CATEGORY_COUNT, MAX_CATEGORY_COUNT );
	_maccum = mix_create( max_sample_count,   MAX_CATEGORY_COUNT );
	_naccum = con_create( max_sample_count );

	_Lwaste = mix_create( max_sample_count,   MAX_CATEGORY_COUNT );
	_Rwaste = mix_create( max_sample_count,   MAX_CATEGORY_COUNT );

	cat_setMinCellCount( _caccum, arg_min_cell_count );

	// DON'T register covan_fini call here because Python extension
	// will do explicit covan_fini; executable's main must register 
	// covan_fini.

	return NULL == _caccum
		|| NULL == _maccum 
		|| NULL == _naccum
		|| NULL == _Lwaste
		|| NULL == _Rwaste ? -1 : 0;
}


/**
 * The input arrays are typed as floats to simplify NaN detection, but this
 * code relies entirely on the equality in sizeof(unsigned int) and 
 * sizeof(float)...and either or both could in fact be const unsigned int*.
 *
 * The following code is slightly longer than strictly necessary because I 
 * want to accomplish everything required to set up subsequent statistical 
 * analyses in ONE ITERATION of the data arrays.
 *
 * There are two major things that happen in each for loop below:
 * 1) Pairs of input covariates are FILTERED.
 *    Only those scalar pairs that are both non-missing are copied
 *    ("pushed") into the class instance that will do the analysis.
 * 2) Arrays are populated for a secondary analysis (Kruskal-Wallis) to
 *    determine whether the missing/non-missing values in one feature
 *    meaningfully segregate the present values in the parallel 
 *    feature (something we would prefer NOT happen).
 *
 * Analyses:
 * 1. both continuous          => Spearman correlation
 * 2. one continuous, one cat. => 1-way ANOVA (or Kruskal-Wallis)
 * 3. both categorical         => contingency table of some kind
 *    3a. both boolean         => Fisher exact
 *    3b. at least one n-ary   => Chi-square
 * 4. Check for statistical difference between each feature segregated
 *    into two groups according to the NA state of its covariate.
 * 5. An "auxiliary" Spearman rho is always computed unless one of 
 *    covariates is categorical with > 2 categories.
 */
int covan_exec( 
		const struct feature_pair *pair,
		struct CovariateAnalysis *covan ) {

	const unsigned int C1
		= pair->l.prop.categories;
	const unsigned int C2
		= pair->r.prop.categories;

	unsigned unused1 = 0;
	unsigned unused2 = 0;
	unsigned count   = 0;

	/**
	  * Insure all string args are initialized to -something- so that 
	  * emitters need be slowed by pervasive NULL checks...
	  */
	covan->result.name
		= covan->waste[0].result.name
		= covan->waste[1].result.name
		= "?";
	// TODO: Revisit this... 
	strcpy( covan->waste[0].result.log, "-" );
	strcpy( covan->waste[1].result.log, "-" );
	strcpy( covan->result.log, "-" );

	/**
	  * Determine feature classes and check for univariate degeneracy.
	  * Univariate degeneracy precludes further analysis. 
	  */

	if( pair->l.prop.constant || pair->r.prop.constant ) {
		covan->status = COVAN_E_UNIVAR_DEGEN;
		return -1;
	} else
	if( C1 > MAX_CATEGORY_COUNT || C2 > MAX_CATEGORY_COUNT ) {
		covan->status = COVAN_E_TOOMANY_CATS;
		return -1;
	}

	// Because we don't anticipate ordinal features yet...

	assert( (pair->l.prop.integral  > 0) == (C1 > 0) );
	assert( (pair->r.prop.integral > 0) == (C2 > 0) );

	covan->stat_class.left 
		= C1 > 0 ? MTM_STATCLASS_CATEGORICAL : MTM_STATCLASS_CONTINUOUS;
	covan->stat_class.right
		= C2 > 0 ? MTM_STATCLASS_CATEGORICAL : MTM_STATCLASS_CONTINUOUS;

	// At this point there should be no other returns until function's end!
	// Collect and report whatever we can...

	mix_clear( _Lwaste, 2 ); // Secondary analyses ALWAYS involve...
	mix_clear( _Rwaste, 2 ); // ...only categories {0,1}.

	// No matter what tests are executed there are only three fundamental
	// cases:
	// 1. categorical-categorical
	// 2. categorical-continuous
	// 3. continuous-continuous
	// Ordinal, which would double the cases, are not currently handled.

	if( covan->stat_class.left == covan->stat_class.right ) {

		if( covan->stat_class.left == MTM_STATCLASS_CONTINUOUS ) {

			con_clear( _naccum );

			for(int i = 0; i < max_sample_count; i++ ) {

				const float F1
					= ((const float*)pair->l.data)[i];
				const float F2
					= ((const float*)pair->r.data)[i];

				if( ! isnan(F1) ) {
					if( ! isnan(F2) ) {
						con_push( _naccum, F1, F2 );
						mix_push( _Lwaste, F1, 1 );
						mix_push( _Rwaste, F2, 1 );
					} else {
						mix_push( _Lwaste, F1, 0 );
						unused1++;
					}
				} else { // F1 is N/A
					if( ! isnan(F2) ) {
						mix_push( _Rwaste, F2, 0 );
						unused2++;
					}
				}
			}

			count = con_size( _naccum );

			if( ! con_complete( _naccum ) ) {
				covan->status |= COVAN_E_COVAR_DEGEN;
			} else
			if( count >= arg_min_sample_count ) {
				con_spearman_correlation( _naccum, &covan->result );
				// TODO: Following line won't be necessary after output formatting is re-implemented for V2.0.
				covan->sign = covan->result.value;
			} else
				covan->status |= COVAN_E_SAMPLES_SIZE;

		} else {

			assert( covan->stat_class.left == MTM_STATCLASS_CATEGORICAL );

			cat_clear( _caccum, C1, C2 );

			// Note that cardinality of univariate features says NOTHING about 
			// the final table after pairs with NA's are removed; it could be
			// empty! That will fall out below though.

			for(int i = 0; i < max_sample_count; i++ ) {

				const unsigned int F1 
					= pair->l.data[i];
				const unsigned int F2 
					= pair->r.data[i];

				// Notice: Though we can't (currently) statistically compare the
				// within-feature discrepancy in CC case as in case involving a
				// continuous variable, I nonetheless use the output structs to
				// at least count the number of wasted samples in each feature.

				if( NAN_AS_UINT != F1 ) {
					if( NAN_AS_UINT != F2 ) {
						cat_push( _caccum, F1, F2 );
					} else { 
						unused1++;
					}
				} else { // F1 is N/A
					if( NAN_AS_UINT != F2 ) {
						unused2++;
					}
				}
			}

			count = cat_size( _caccum );

			if( ! cat_complete( _caccum ) ) {
				covan->status |= COVAN_E_COVAR_DEGEN;
			} else {
				cat_cullBadCells( _caccum, covan->result.log, MAXLEN_STATRESULT_LOG );
				// ...cullBadCells won't allow the table to become degenerate. 
				if( count >= arg_min_sample_count ) {
					if( cat_is2x2( _caccum ) ) {
						cat_fisher_exact( _caccum, &covan->result );
						covan->sign = cat_spearman_rho( _caccum );
					} else {
						cat_chi_square( _caccum, &covan->result );
					}
				} else
					covan->status |= COVAN_E_SAMPLES_SIZE;
			}
		}

	} else { // features are not of same class

		if( covan->stat_class.left == MTM_STATCLASS_CATEGORICAL ) {

			assert( covan->stat_class.right == MTM_STATCLASS_CONTINUOUS );

			mix_clear( _maccum, C1 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const unsigned int F1 
					= pair->l.data[i];
				const float F2 
					= ((const float*)pair->r.data)[i];

				if( ! isnan(F2) ) {
					if( NAN_AS_UINT != F1 ) {
						mix_push( _maccum, F2, F1 );
						mix_push( _Rwaste, F2, 1 );
					} else {
						mix_push( _Rwaste, F2, 0 );
						unused2++;
					}
				} else { // F2 is N/A
					if( NAN_AS_UINT != F1 ) {
						unused1++;
					}
				}
			}

		} else {

			assert( covan->stat_class.left == MTM_STATCLASS_CONTINUOUS && 
					covan->stat_class.right == MTM_STATCLASS_CATEGORICAL );

			mix_clear( _maccum, C2 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const float F1 
					= ((const unsigned int*)pair->l.data)[i];
				const unsigned int F2 
					= pair->r.data[i];

				if( ! isnan(F1) ) {
					if( NAN_AS_UINT != F2 ) {
						mix_push( _maccum, F1, F2 );
						mix_push( _Lwaste, F1, 1 );
					} else {
						mix_push( _Lwaste, F1, 0 );
						unused1++;
					}
				} else { // F1 is N/A
					if( NAN_AS_UINT != F2 ) {
						unused2++;
					}
				}
			}
		}

		count = mix_size( _maccum );

		if( ! mix_complete( _maccum ) ) {
			covan->status |= COVAN_E_COVAR_DEGEN;
		} else
		if( count >= arg_min_sample_count ) {
			mix_kruskal_wallis( _maccum, &covan->result );
			if( mix_categoricalIsBinary( _maccum ) )
				covan->sign = mix_spearman_rho( _maccum );
		} else
			covan->status |= COVAN_E_SAMPLES_SIZE;
	}

	/**
	  * Some tests are ordered
	  */

	/**
	  * post-filtering counts should always be valid whatever 
	  * else has occurred.
	  */
	covan->waste[0].unused = unused1;
	covan->waste[1].unused = unused2;

	// Characterize how the unused parts of the two samples might have
	// affected the statistics computed on their "overlap".

	if( mix_complete( _Lwaste ) )
		mix_kruskal_wallis( _Lwaste, &(covan->waste[0].result) );

	if( mix_complete( _Rwaste ) )
		mix_kruskal_wallis( _Rwaste, &(covan->waste[1].result) );

	return covan->status ? -1 : 0;
}

