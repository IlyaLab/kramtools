
/**
  *
  * For reference: Pearson's rho...
  *
  *               sum( (x_i-E[x])(y_i-E[y]) )
  * rho = -------------------------------------------
  *       sqrt( sum((x_i-E[x])^2) sum((y_i-E[y])^2) )
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <alloca.h>

#include <gsl/gsl_cdf.h>

#include "stattest.h"
#include "mix.h"
#include "bvr.h"
#include "limits.h"

struct Pair {
	float cv; // Continuous Variable
	unsigned int dv; // Discrete Variable
};


static int _cmp_pair( const void *pvl, const void *pvr ) {
	const struct Pair *l = (const struct Pair *)pvl;
	const struct Pair *r = (const struct Pair *)pvr;
	const float delta = l->cv - r->cv;
	if( delta == 0.0 )
		return 0.0;
	else
		return  delta < 0.0 ? -1 : +1;
}


struct MixCovars {

	/**
	  * Maximum number of samples
	  */
	unsigned int SAMPLE_CAPACITY;

	/**
	  * All category labels pushed into this accumulator must
	  * be in [0,category_capacity).
	 */ 
	unsigned int CATEGORY_CAPACITY;

	/**
	  * Count of samples pushed since last cat_clear.
	  */
	unsigned int sample_count;

	/**
	  * This essentially bounds the allowed category LABELS.
	  * Any category label push'ed in must be in [0,expected_categories).
	  * Only this member, NOT category_capacity, should be used to size 
	  * counts arrays!
	  */
	unsigned int expected_categories;

	/**
	  * This is an count of DISTINCT LABELS observed by push.
	  */
	unsigned int observed_categories;

	// Calculation state ///////////////////////////////////////////////////

	double mean_rank;  // necessarily of both covariates
	double sum_sq_dev; // used for Kruskal-Wallis and Spearman

	/**
	  * Following struct array is ONLY for Spearman rho calculation.
	  * Note that the number of non-empty categories after data entry
	  * says NOTHING about WHICH categories (indices) are non-zero.
	  * An earlier implementation was assuming category_count[0,1] > 0,
	  * BUT THIS NEED NOT BE TRUE! Thus, all this hack...
	  */
	struct {
		unsigned int index;
		unsigned int count;
		unsigned int meanRank;
	} edge[2];

	double sum_dev_prod; // used only for Spearman-rho

	// Buffer //////////////////////////////////////////////////////////////

	/**
	  * The TOTAL amount of malloc'ed space (for the purposes
	  * of fast clearing with memset).
	  */
	size_t SIZEOF_BUFFERS;

	/**
	  * Two buffers
	  */
	struct Pair *samples;

	/**
	  * A buffer allocated once in construction and reused for 
	  * all subsequent analyses.
	  * Reallocation is not currently supported.
	  */
	unsigned int *category_count;
};


#ifdef _UNITTEST_MIX_
unsigned arg_min_mixb_count = 1;
#else
extern unsigned arg_min_mixb_count;
#endif

/**
  * Both minCat and maxCat assume that samples is non-empty,
  * but it's entirely possible that there were NO samples.
  */
static unsigned _minCat( struct MixCovars *co ) {
	unsigned int i = 0;
	while( i < co->expected_categories ) {
		if( co->category_count[i] > 0 ) 
			return i;
		i++;
	}
	return INT_MAX; // just to silence compiler.
}


static unsigned _maxCat( struct MixCovars *co ) {

	int i = co->expected_categories;
	while( i-- > 0 ) {
		if( co->category_count[i] > 0 ) 
			return (unsigned int)i;
	}
	return INT_MAX; // just to silence compiler
}


/**
  * Notice that all calcuations are performed without actually requiring
  * an array of ranks.
  */
static unsigned int _rank_sums( struct MixCovars *co,
	double *sums ) {

	const unsigned int N = co->sample_count;
	const bool PRECALC_SPEARMAN_PRODUCT
		= co->observed_categories == 2;
	const unsigned UPPER
		= co->edge[1].index;
	unsigned int ties = 0;
	unsigned int until = 0;
	float rank = 0;
	double diff;

	// The default pair::operator< automatically uses the first element.

	qsort( co->samples, co->sample_count, sizeof(struct Pair), _cmp_pair );

	for(unsigned int i = 0; i < N; i++ ) {

		const unsigned int cat 
			= co->samples[i].dv;

		// The rank of the current sample is i+1 UNLESS we're in
		// the midst of a run of samples with the same value in
		// which case the last computed rank (a mean) is reused.

		if( ! ( i < until ) ) {

			until = i + 1;
			while( until < N 
					&& ( co->samples[i].cv == co->samples[until].cv ) ) {
				ties++;
				until++;
			}

			if( i == 0 && until == N ) {
				// degenerate case: all ties.
			}

			// Mean of n integers starting on r:
			// [ (r+0) + (r+1) + (r+2) + ... + (r+(n-1)) ] / n
			// ... = \frac{nr + \sum_{i=0}^{n-1} i}{n} 
			// ... = r + (1/n)[ (n-1)n / 2 ]
			// ... = r + (n-1)/2

			rank = (1.0+i) + ((until-i)-1)/2.0;

			// ...notice if the very next sample -is- different
			// then rank == i+1, the trivial case.
		}

		/**
		  * Warning: If the numeric feature is constant then rank==meanRank
		  * and diff == 0
		  */
		diff = rank - co->mean_rank;

		co->sum_sq_dev += ( diff * diff );

		if( PRECALC_SPEARMAN_PRODUCT ) {
			co->sum_dev_prod += ( diff * ( co->edge[UPPER==cat?1:0].meanRank - co->mean_rank ) );
		}

		sums[ cat ] += rank;
	}

	return ties;
}


#if defined(_UNITTEST_MIX_)
static void dbg_dump( struct MixCovars *co, FILE *fp ) {
	for(unsigned int i = 0; i < co->sample_count; i++ )
		fprintf( fp, "%.3f\t%d\n", co->samples[i].cv, co->samples[i].dv );
}
#endif

/***************************************************************************
  * Publics
  */

void mix_destroy( void *pv ) {

	if( pv ) {
		struct MixCovars *co = (struct MixCovars *)pv;
		if( co->samples )
			free( co->samples );
		free( pv );
	}
}


/**
  * Pre-allocate a set of working buffers large enough for all anticipated
  * calculations (max feature length) and a struct to wrap them.
  */
void *mix_create( unsigned int sample_capacity, unsigned int category_capacity ) {

	struct MixCovars *co
		= calloc( 1, sizeof(struct MixCovars) );
	if( co ) {
		co->SAMPLE_CAPACITY   = sample_capacity;
		co->CATEGORY_CAPACITY = category_capacity;
		co->SIZEOF_BUFFERS
			= ( sample_capacity * sizeof(struct Pair) )
			+ ( category_capacity * sizeof(unsigned int) );

		// Allocate one large buffer and partition it up.
		co->samples
			= calloc( co->SIZEOF_BUFFERS, sizeof(char) );
		co->category_count
			= (unsigned int*)(co->samples
			+ sample_capacity);
		// If -anything- failed clean up any successes.
		if( NULL == co->samples ) {
			mix_destroy( co );
			return NULL;
		}
		return co;
	}
	return NULL;
}


void mix_clear( void *pv, unsigned expcat ) {

	struct MixCovars *co = (struct MixCovars *)pv;

	assert( expcat <= co->CATEGORY_CAPACITY );

	co->sample_count = 0;
	co->expected_categories = expcat;
	co->observed_categories = 0;
	co->mean_rank    = 0.0;
	co->sum_sq_dev   = 0.0;
	co->sum_dev_prod = 0.0;
	memset( co->edge, 0, sizeof(co->edge) );
	co->sum_dev_prod = 0.0;

	memset( co->samples, 0, co->SIZEOF_BUFFERS );
}


void mix_push( void *pv, float num, unsigned int cat ) {

	struct MixCovars *co = (struct MixCovars *)pv;

	assert( cat              < co->expected_categories );
	assert( co->sample_count < co->SAMPLE_CAPACITY     );

	if( 0 == co->category_count[ cat ] ) 
		co->observed_categories++;

	co->category_count[ cat ] += 1;

	// Not updating edges in here because the number of conditionals
	// executed for sample counts > 32 exceeds the work to find the
	// edges post-sample accumulation.
	co->samples[ co->sample_count ].cv = num;
	co->samples[ co->sample_count ].dv = cat;
	co->sample_count++;
}


size_t mix_size( void *pv ) {
	struct MixCovars *co = (struct MixCovars *)pv;
	return co->sample_count;
}


/**
  * - Establish whether or not Spearman is even sensible (according to
  *   whether or not the categorical variable is binary).
  * - Compute the mean rank used throughout subsequent calcs.
  * - Compute the mean ranks of the categorical values assuming
  *   they -are- binary.
  */
bool mix_complete( void *pv ) {

	struct MixCovars *co = (struct MixCovars *)pv;
	if( mix_degenerate( co ) ) return false;

	co->mean_rank = (1.0 + co->sample_count)/2.0;

	co->edge[1].index = _maxCat( co ); // ...this is used by kruskal_wallis, but

	assert( co->edge[1].index < co->expected_categories );

	// ...the remainder of the edge[] struct is only used if...

	if( co->observed_categories == 2 ) {

		// All the following is relevant only to Spearman rho calculation.

		co->edge[0].index = _minCat( co );
		assert( co->edge[0].index < co->expected_categories );

		co->edge[0].count = co->category_count[ co->edge[0].index ];
		co->edge[0].meanRank = MEAN_RANK_OF_TIES( 0, co->edge[0].count );

		co->edge[1].count = co->category_count[ co->edge[1].index ];

		if( co->edge[0].count < arg_min_mixb_count ||
			co->edge[1].count < arg_min_mixb_count )
			return false; // another degeneracy class

		co->edge[1].meanRank = MEAN_RANK_OF_TIES( co->edge[0].count, co->sample_count );

		assert( (co->edge[0].count + co->edge[1].count) == co->sample_count );
	}
	return true;
}


bool mix_degenerate( void *pv ) {
	struct MixCovars *co = (struct MixCovars *)pv;
	return co->observed_categories < 2 || co->sample_count < 2;
}


bool mix_categoricalIsBinary( void *pv ) {
	struct MixCovars *co = (struct MixCovars *)pv;
	return co->observed_categories == 2;
}


/**
 * Calculates the Kruskal-Wallis statistic and optionally a p-value.
 * Importantly, it does it in a one pass iteration over the data.
 */
int mix_kruskal_wallis( void *pv, struct Statistic *result ) {

	struct MixCovars *co = (struct MixCovars *)pv;
	const unsigned int N 
		= co->sample_count ;
	const size_t SIZEOF_SUMS
		= co->expected_categories * sizeof(double);
	double *rank_sum
		= (double*)alloca( SIZEOF_SUMS );
	double numerator = 0.0;

	memset( rank_sum, 0, SIZEOF_SUMS   );

	result->extra_value[0] = _rank_sums( co, rank_sum );

	for(unsigned int i = 0; i <= co->edge[1].index; i++ ) {
		if( co->category_count[i] > 0 ) {
			rank_sum[i] /= co->category_count[i];
			double delta 
				= rank_sum[i] - co->mean_rank;
			numerator 
				+= ( co->category_count[i] * delta * delta );
		}
	}

	result->name
		= "Kruskal-Wallis_K";
	result->sample_count
		= N;
	result->value
		= (N-1) * ( numerator / co->sum_sq_dev );
	result->probability
		= gsl_cdf_chisq_Q( result->value, co->observed_categories-1 );

	return 0;
}


/**
 * P-value returned is ONE-SIDED. Double it for two-sided.
 */
#ifdef HAVE_MANN_WHITNEY
int mix_mann_whitney( void *pv, struct Statistic *result ) {

	struct MixCovars *co = (struct MixCovars *)pv;
	const size_t SIZEOF_SUMS
		= co->expected_categories * sizeof(double);
	double *rank_sum
		= (double*)alloca( SIZEOF_SUMS );
	memset( rank_sum, 0, SIZEOF_SUMS   );

	ms->ties
		= _rank_sums( co, rank_sum );

	// TODO: If too many ties, the following method of calculating U
	//       becomes invalid. What to do?

	const unsigned int SQ
		= category_count[0]*category_count[1];

	// If either group count is 0 we're done...

	if( SQ ) {

		const double U1
			= rank_sum[0] - (category_count[0]*(category_count[0]+1.0))/2.0;
		const double U2
			= rank_sum[1] - (category_count[1]*(category_count[1]+1.0))/2.0;

		assert( U1+U2 == SQ );

		ms->U = U1 < U2 ? U1 : U2;

		const double mean_U
			= SQ / 2.0;
		const double sigma_U
			= sqrt( SQ*( category_count[0] + category_count[1] + 1.0 ) / 12.0 );

		cs->P = gsl_cdf_ugaussian_P( ( ms->U - mean_U ) / sigma_U );
	}
	return 0;
}
#endif

/**
  * This calculation depends on a quantities calculated within
  * kruskal_wallis (or mann_whitney), so must be called after
  * either of them.
  * (This is also the only reason this method can be const.)
  */
double mix_spearman_rho( void *pv ) {

	struct MixCovars *co = (struct MixCovars *)pv;
	/**
	  * denominator requires sum of square deviation for
	  * the binary variables which collapses to...
	  */
	double SQD0 = (co->edge[0].meanRank - co->mean_rank ); SQD0 *= SQD0;
	double SQD1 = (co->edge[1].meanRank - co->mean_rank ); SQD1 *= SQD1;

	const double SUMSQD
		= co->edge[0].count * SQD0 + co->edge[1].count * SQD1;
	/**
	  * And the sum of squared deviation of the numeric covariate's
	  * ranks (sum_sq_dev) was (necessarily) calculated earlier!
	  */
	return co->sum_dev_prod / sqrt( co->sum_sq_dev * SUMSQD );
}


#ifdef _UNITTEST_MIX_

/**
 * This is intended to be exercised with the following R script
 * that generates a small table of grouped floating-point values
 * executes a Kruskal-Wallis test on it, and dumps the table
 * to a tab-delimited file with the test results as a comment on
 * the first line.
 *
 * x <- data.frame(
 * 		cat=as.integer(gl(3,4))-1, 
 * 		num=c( runif(4)*10, runif(4)*20, 10+runif(4)*10 ) );
 * k <- with( x, kruskal.test( num, cat ) );
 * cat( sprintf( "# K=%f p-value=%f\n", k$statistic, k$p.value ), file="foo.tab" );
 * write.table( x, 'foo.tab', quote=F, sep='\t', row.names=F, col.names=F, append=TRUE );
 */
#include <err.h>

int main( int argc, char *argv[] ) {

	if( argc >= 3 ) {

		struct Statistic result;
		const int EXPCAT
			= atoi( argv[1] );
		struct MixCovars *accum
		   = mix_create( atoi( argv[2] ), MAX_CATEGORY_COUNT );
		FILE *fp 
			= argc > 3
			? fopen( argv[3], "r" )
			: stdin;
		char *line = NULL;
		size_t n = 0;

		mix_clear( accum, EXPCAT );

		while( getline( &line, &n, fp ) > 0 ) {
			unsigned int cat;
			float        num;
			if( line[0] == '#' ) {
				fputs( line, stdout );
				continue;
			}
			if( 2 == sscanf( line, "%d\t%f\n", &cat, &num ) ) {
				mix_push( accum, num, cat );
			} else {
				fprintf( stderr, "Failure parsing line: %s", line );
			}
		}
		free( line );
		fclose( fp );

		// Always calculate Kruskal-Wallis since it's valid for any groups
		// count >= 2...

		if( mix_complete( accum ) ) {
			memset( &result, 0, sizeof(struct Statistic) );
			if( mix_kruskal_wallis( accum, &result ) == 0 ) {
				printf( "K=%f, p-value=%f", result.value, result.probability );
				if( mix_categoricalIsBinary( accum ) )
					printf( ", spearman=%f\n", mix_spearman_rho(accum) );
				else
					fputc( '\n', stdout );
			} else {
				printf( "error\n" );
			}
		}

#ifdef HAVE_MANN_WHITNEY
		// Do a Mann-Whitney, too, if there are only 2 groups.
		if( mix_categoricalIsBinary(accum) == 2 ) {
		}
#endif
		mix_destroy( accum );
	} else
		err( -1, "%s <categories> <sample count> [ <input file> ]", argv[0] );
	return 0;
}
#endif

