/**
 * Three C++ classes handle all actual analysis:
 * 1) CatCovars deals with pairs of categorical vectors.
 * 2) MixCovars deals with pairs containing one numeric and one categoric
 *    variable.
 * 3) NumCovars deals with pairs of numeric vectors.
 *
 * This module is really just a dispatcher and metadata collector for
 * all four possibilities:
 * 1) numeric, numeric
 * 2) categorical, numeric
 * 3) numeric, categorical 
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

#include <ostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>

#include "stattest.h"
#include "analysis.h"
#include "cat.h"
#include "mix.h"
#include "num.h"
#include "args.h"
#include "mtmatrix.h"

/**
  * Most arrays are pre-allocated and sized according to this variable.
  */
static int max_sample_count = 0;

/**
 * These classes handle the actual feature1 vs feature2 analyses.
 */
static CatCovars _caccum( _MAX_CATEGORIES, _MAX_CATEGORIES );
static MixCovars _maccum( _MAX_CATEGORIES );
static NumCovars _naccum;

/**
 * These classes handle comparisons between the two groups WITHIN a
 * single feature distinguished according to whether the corresponding
 * covariate (in the other feature) is present (not "NA").
 * In numeric/numeric comparisons are both used, in numeric/categorical
 * comparisons one or the other is used, and in cat/cat neither is used.
 *
 * Notice that in general counts by category in categorical/numeric
 * pairs are generally unavailable with one exception: the category 0's
 * count is always output by the Kruskal-Wallis test. This is motivated
 * primarily by the application defined above. As a result, it is 
 * ESSENTIAL THAT SAMPLES PRESENT IN FEAT1 BUT MISSING IN FEAT2 ARE
 * TAGGED WITH 0 IN _Lwaste (vica versa for 2).
 */
static MixCovars _Lwaste( _MAX_CATEGORIES );
static MixCovars _Rwaste( _MAX_CATEGORIES );

/**
 * Let C++ choose the appropriate one...
 */
inline bool _present( unsigned int u ) {
	return 0x7FC00000 != u;
}

inline bool _present( float f ) {
	return not isnan(f);
}

////////////////////////////////////////////////////////////////////////////
// Public API
////////////////////////////////////////////////////////////////////////////

/**
 * This must pre-allocate all the memory we might need for every combination
 * of data types.
 */
int covan_init( int columns ) {

	max_sample_count = columns; // EVERYTHING depends on this.

	_naccum.reserve( max_sample_count );
	_maccum.reserve( max_sample_count );
	_Lwaste.reserve( max_sample_count );
	_Rwaste.reserve( max_sample_count );

	_caccum.setMinCellCount( arg_min_cell_count );

	return 0;
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
 * 1. both numeric            => Spearman correlation
 * 2. one numeric, one cat.   => 1-way ANOVA (or Kruskal-Wallis)
 * 3. both categorical        => contingency table of some kind
 *    3a. both boolean        => Fisher exact
 *    3b. at least one n-ary  => Chi-square
 * 4. Check for statistical difference between each feature segregated
 *    into two groups according to the NA state of its covariate.
 * 5. An "auxiliary" Spearman rho is always computed unless one of 
 *    covariates is categorical with > 2 categories.
 */
int covan_exec( 
		const struct mt_row_pair *pair,
		struct CovariateAnalysis *covan ) {

	const unsigned int C1
		= pair->l.prop.categories;
	const unsigned int C2
		= pair->r.prop.categories;

	unsigned unused1 = 0;
	unsigned unused2 = 0;
	unsigned count   = 0;

	/**
	  * Determine feature classes and check for univariate degeneracy.
	  * Univariate degeneracy precludes further analysis. 
	  */

	if( pair->l.prop.constant || pair->r.prop.constant ) {
		covan->status = FAIL_E_DEGEN;
		return -1;
	} else
	if( C1 > _MAX_CATEGORIES || C2 > _MAX_CATEGORIES ) {
		covan->status = FAIL_TOOMANY;
		return -1;
	}

	// Because we don't anticipate ordinal features yet...

	assert( (pair->l.prop.integral  > 0) == (C1 > 0) );
	assert( (pair->r.prop.integral > 0) == (C2 > 0) );

	covan->lclass = C1 > 0 ? Categorical : Continuous;
	covan->rclass = C2 > 0 ? Categorical : Continuous;

	// At this point there should be no other returns until function's end!
	// Collect and report whatever we can...

	_Lwaste.clear( 2 ); // Secondary analyses ALWAYS involve...
	_Rwaste.clear( 2 ); // ...only categories {0,1}.

	// No matter what tests are executed there are only three fundamental
	// cases:
	// 1. categorical-categorical
	// 2. categorical-continuous
	// 3. continuous-continuous
	// Ordinal, which would double the cases, are not currently handled.

	if( covan->lclass == covan->rclass ) {

		if( covan->lclass == Continuous ) {

			_naccum.clear();

			for(int i = 0; i < max_sample_count; i++ ) {

				const float F1
					= reinterpret_cast<const float*>(pair->l.data)[i];
				const float F2
					= reinterpret_cast<const float*>(pair->r.data)[i];

				if( _present(F1) ) {
					if( _present(F2) ) {
						_naccum.push( F1, F2 );
						_Lwaste.push( F1, 1 );
						_Rwaste.push( F2, 1 );
					} else {
						_Lwaste.push( F1, 0 );
						unused1++;
					}
				} else { // F1 is N/A
					if( _present(F2) ) {
						_Rwaste.push( F2, 0 );
						unused2++;
					}
				}
			}

			count = _naccum.size();

			if( not _naccum.complete() ) {
				covan->status |= FAIL_L_DEGEN;
			} else
			if( count >= arg_min_sample_count ) {
				_naccum.spearman_correlation( &covan->result );
			} else
				covan->status |= FAIL_SAMPLES;

		} else {

			assert( covan->lclass == Categorical );

			_caccum.clear( C1, C2 );

			// Note that cardinality of univariate features says NOTHING about 
			// the final table after pairs with NA's are removed! It could be
			// empty! That will fall out below though.

			for(int i = 0; i < max_sample_count; i++ ) {

				const unsigned int F1 
					= pair->l.data[i];
				const unsigned int F2 
					= pair->r.data[i];

				// Notice: Though we can't (currently) statistically compare the
				// within-feature discrepancy in CC case as in case involving a
				// numeric variable, I nonetheless use the output structs to
				// at least count the number of wasted samples in each feature.

				if( _present(F1) ) {
					if( _present(F2) ) {
						_caccum.push( F1, F2 );
					} else { 
						unused1++;
					}
				} else { // F1 is N/A
					if( _present(F2) ) {
						unused2++;
					}
				}
			}

			count = _caccum.size();

			if( not _caccum.complete() ) {
				covan->status |= FAIL_L_DEGEN;
			} else {
				_caccum.cullBadCells( covan->result.log, MAXLEN_STATRESULT_LOG );
				if( count >= arg_min_sample_count ) {
					if( _caccum.is2x2() ) {
						_caccum.fisher_exact( &covan->result );
						covan->sign = 0;// TODO: sign()
					} else {
						_caccum.chi_square( &covan->result );
					}
				} else
					covan->status |= FAIL_SAMPLES;
			}
		}

	} else { // features are not of same class

		if( covan->lclass == Categorical ) {

			assert( covan->rclass == Continuous );

			_maccum.clear( C1 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const unsigned int F1 
					= pair->l.data[i];
				const float F2 
					= reinterpret_cast<const float*>(pair->r.data)[i];

				if( _present(F2) ) {
					if( _present(F1) ) {
						_maccum.push( F2, F1 );
						_Rwaste.push( F2, 1 );
					} else {
						_Rwaste.push( F2, 0 );
						unused2++;
					}
				} else { // F2 is N/A
					if( _present(F1) ) {
						unused1++;
					}
				}
			}

		} else {

			assert( covan->lclass == Continuous && covan->rclass == Categorical );

			_maccum.clear( C2 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const float F1 
					= reinterpret_cast<const unsigned int*>(pair->l.data)[i];
				const unsigned int F2 
					= pair->r.data[i];

				if( _present(F1) ) {
					if( _present(F2) ) {
						_maccum.push( F1, F2 );
						_Lwaste.push( F1, 1 );
					} else {
						_Lwaste.push( F1, 0 );
						unused1++;
					}
				} else { // F1 is N/A
					if( _present(F2) ) {
						unused2++;
					}
				}
			}
		}

		count = _maccum.size();

		if( not _maccum.complete() ) {
			covan->status |= FAIL_L_DEGEN;
		} else
		if( count >= arg_min_sample_count ) {
			_maccum.kruskal_wallis( &covan->result );
			if( _maccum.categoricalIsBinary() )
				covan->sign = 0;// TODO: sign()
		} else
			covan->status |= FAIL_SAMPLES;
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

	if( _Lwaste.complete() )
		_Lwaste.kruskal_wallis( &(covan->waste[0].result) );

	if( _Rwaste.complete() )
		_Rwaste.kruskal_wallis( &(covan->waste[1].result) );

	return covan->status ? -1 : 0;
}

