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

/**
 * These constants must match similar constants used in the Python
 * script that converts TCGA tabbed-separated-values files into the
 * binary form that this code expects. (...currently named prep.py)
 */
#define MISSING_MASK     0x0000FFFF
#define CATCOUNT_MASK    0x00FF0000
#define CATEGORICAL_FLAG 0x01000000

/**
  * Most arrays are pre-allocated and sized according to this variable.
  */
static int max_sample_count = 0;

/**
 * These classes handle the actual feature1 vs feature2 analyses.
 */
static CatCovars _caccum( _MAX_CATEGORIES, _MAX_CATEGORIES );
static MixCovars _maccum( _MAX_CATEGORIES, arg_min_mixb_count );
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
 * TAGGED WITH 0 IN _waste1 (vica versa for 2).
 */
static MixCovars _waste1( _MAX_CATEGORIES, arg_min_mixb_count );
static MixCovars _waste2( _MAX_CATEGORIES, arg_min_mixb_count );

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
int analysis_init( int columns ) {

	max_sample_count = columns; // EVERYTHING depends on this.

	_naccum.reserve( max_sample_count );
	_maccum.reserve( max_sample_count );
	_waste1.reserve( max_sample_count );
	_waste2.reserve( max_sample_count );

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
unsigned analysis_exec( 
		unsigned int h1, const float *f1,  // feature 1
		unsigned int h2, const float *f2,  // feature 2
		struct CovarsSummary *res ) {

	static const char *TOO_MANY_CATS
		= "error: category count (%d) in feature %d exceeds maximum (%d)\n";

	// Extract the category counts from the (32-bit unsigned int) header.
	// Bitfields are defined in the Python3 code prep.py.

	const unsigned int C1 = (h1 & CATCOUNT_MASK) >> 16;
	const unsigned int C2 = (h2 & CATCOUNT_MASK) >> 16;

	unsigned status  = 0;
	unsigned unused1 = 0;
	unsigned unused2 = 0;
	unsigned count   = 0;

	/**
	 * Regardless of the classes of variables, if either C1 or C2 equals 
	 * one, we're in a degenerate situation. Either:
	 * 1) a constant numeric vector (i.e. one with a single non-missing
	 *    repeated value), or
	 * 2) a categorical vector with a single (non-missing) category.
	 * With univariate degeneracy there's no point in going further...
	 */
	if( C1 == 1 ) {
		res->kind = DegenIndet;
		return FAIL_E_DEGEN;
	} else
	if( C2 == 1 ) {
		res->kind = IndetDegen;
		return FAIL_E_DEGEN;
	} else 
	if( not ( C1 <= _MAX_CATEGORIES ) ) {
		fprintf( stderr, TOO_MANY_CATS, C1, 1, _MAX_CATEGORIES );
		return FAIL_TOOMANY;
	} else
	if( not ( C2 <= _MAX_CATEGORIES ) ) {
		fprintf( stderr, TOO_MANY_CATS, C2, 2, _MAX_CATEGORIES );
		return FAIL_TOOMANY;
	} else {
		res->kind = (enum CovarTypes)(
		     ((h1 & CATEGORICAL_FLAG)?1:0)
		   + ((h2 & CATEGORICAL_FLAG)?2:0) 
		   + 1 );
		assert( UnknownCovar < res->kind and res->kind <= CatCat );
	}

	// At this point there should be no other returns until function's end!
	// Collect and report whatever we can...

	_waste1.clear( 2 ); // Secondary analyses ALWAYS involve...
	_waste2.clear( 2 ); // ...only categories {0,1}.

	switch( res->kind ) {

	case NumNum:

		_naccum.clear();

		for(int i = 0; i < max_sample_count; i++ ) {

			const float F1 = f1[i];
			const float F2 = f2[i];

			if( _present(F1) ) {
				if( _present(F2) ) {
					_naccum.push( F1, F2 );
					_waste1.push( F1, 1 );
					_waste2.push( F2, 1 );
				} else {
					_waste1.push( F1, 0 );
					unused1++;
				}
			} else { // F1 is N/A
				if( _present(F2) ) {
					_waste2.push( F2, 0 );
					unused2++;
				}
			}
		}

		count = _naccum.size();

		if( not _naccum.complete() ) {
			status |= FAIL_L_DEGEN;
		} else
		if( count >= arg_min_sample_count ) {
			res->test = Spearman;
			_naccum.spearman_correlation( &res->common, &res->spearman );
			res->spearman_rho = res->spearman.rho;
		} else
			status |= FAIL_SAMPLES;
		break; // end-of-case NumNum

	case CatNum: // F1 categorical, F2 numeric
	case NumCat: // F1 numeric, F2 categorical

		if( CatNum == res->kind ) {

			_maccum.clear( C1 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const unsigned int F1 
					= reinterpret_cast<const unsigned int*>(f1)[i];
				const float F2 = f2[i];

				if( _present(F2) ) {
					if( _present(F1) ) {
						_maccum.push( F2, F1 );
						_waste2.push( F2, 1 );
					} else {
						_waste2.push( F2, 0 );
						unused2++;
					}
				} else { // F2 is N/A
					if( _present(F1) ) {
						unused1++;
					}
				}
			}

		} else { // NumCat

			_maccum.clear( C2 );

			for(int i = 0; i < max_sample_count; i++ ) {

				const float F1 = f1[i];
				const unsigned int F2 
					= reinterpret_cast<const unsigned int*>(f2)[i];

				if( _present(F1) ) {
					if( _present(F2) ) {
						_maccum.push( F1, F2 );
						_waste1.push( F1, 1 );
					} else {
						_waste1.push( F1, 0 );
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
			status |= FAIL_L_DEGEN;
		} else
		if( count >= arg_min_sample_count ) {
			res->test = KruskalWallis;
			_maccum.kruskal_wallis( &res->common, &res->kruskal );
			if( _maccum.categoricalIsBinary() )
				res->spearman_rho = _maccum.spearman_rho();
		} else
			status |= FAIL_SAMPLES;
		break; // end-of-case CatNum, NumCat

	case CatCat: // both categorical

		_caccum.clear( C1, C2 );

		// Note that cardinality of univariate features says NOTHING about 
		// the final table after pairs with NA's are removed! It could be
		// empty! That will fall out below though.

		for(int i = 0; i < max_sample_count; i++ ) {

			const unsigned int F1 
				= reinterpret_cast<const unsigned int*>(f1)[i];
			const unsigned int F2 
				= reinterpret_cast<const unsigned int*>(f2)[i];

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
			status |= FAIL_L_DEGEN;
		} else {
			_caccum.cullBadCells( res->log, SUMMARY_LOG_LEN );
			if( count >= arg_min_sample_count ) {
				if( _caccum.is2x2() ) {
					res->test = FisherExact;
					_caccum.fisher_exact( &res->common, &res->fisherx );
					res->spearman_rho = _caccum.spearman_rho();
				} else {
					res->test = ChiSquare;
					_caccum.chi_square( &res->common, &res->chisq );
				}
			} else
				status |= FAIL_SAMPLES;
		}
		break; // end-of-case CatCat

	default:
		fprintf( stderr, "Should be impossible to reach %s:%d!", 
			__FILE__, __LINE__ );
		abort();
	}

	/**
	  * post-filtering counts should always be valid whatever 
	  * else has occurred.
	  */
	res->common.N        = count;
	res->waste[0].unused = unused1;
	res->waste[1].unused = unused2;

	// Characterize how the unused parts of the two samples might have
	// affected the statistics computed on their "overlap".

	if( _waste1.complete() )
		_waste1.kruskal_wallis( 
			&(res->waste[0].common), 
			&(res->waste[0].kruskal) );

	if( _waste2.complete() )
		_waste2.kruskal_wallis( 
			&(res->waste[1].common), 
			&(res->waste[1].kruskal) );

	return status;
}

