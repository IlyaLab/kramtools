
/**
  *
  * For reference: Pearson's rho...
  *
  *               sum( (x_i-E[x])(y_i-E[y]) )
  * rho = -------------------------------------------
  *       sqrt( sum((x_i-E[x])^2) sum((y_i-E[y])^2) )
  */

#include <ostream>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <alloca.h>
#include <gsl/gsl_cdf.h>

#include "stattest.h"
#include "mix.h"
#include "bvr.h"
#include "limits.h"

#ifdef _UNITTEST_MIX_
unsigned arg_min_mixb_count = 1;
void panic( const char *src, int line ) {
	fprintf( stderr, "panic on %s:%d", src, line );
	abort();
}
#else
extern "C" unsigned arg_min_mixb_count;
extern "C" void panic( const char *src, int line );
#endif

using namespace std;

#if defined(_DEBUG) || defined(_UNITTEST_MIX_)
ostream& operator<<( ostream& os, const MixCovars& ob ) {

	for(unsigned int i = 0; i < ob.size(); i++ ) {
		os << ob.samples[i].first << "\t" << ob.samples[i].second << endl;
	}
	return os;
}
#endif

MixCovars::MixCovars( unsigned int cap ) : 
	CATEGORY_CAPACITY(cap) {
	category_count = new unsigned int[ CATEGORY_CAPACITY ];
	clear(0);
}


MixCovars::~MixCovars() {
	delete [] category_count;
}


void MixCovars::clear( unsigned expcat ) {

	expected_categories = expcat;
	observed_categories = 0;

	assert( expected_categories <= CATEGORY_CAPACITY );

	samples.clear();
	memset( category_count, 0, expected_categories*sizeof(unsigned int)); // OPTIMIZATION
	meanRank         = 0.;
	sum_sq_dev       = 0.;
	// Spearman rho stuff...
	memset( edge, 0, sizeof(edge) ); // confirmed sizeof(edge)==24
	sum_dev_prod     = 0.;
}


void MixCovars::reserve( unsigned int n ) {
	samples.reserve( n ); // just pre-allocates storage.
}


/**
  * Both minCat and maxCat assume that samples is non-empty,
  * but it's entirely possible that there were NO samples.
  */
unsigned MixCovars::minCat() const {
	unsigned int i = 0;
	while( i < expected_categories ) {
		if( category_count[i] > 0 ) 
			return i;
		i++;
	}
	panic( __FILE__, __LINE__ );
	return INT_MAX; // just to silence compiler.
}


unsigned MixCovars::maxCat() const {

	int i = expected_categories;
	while( i-- > 0 ) {
		if( category_count[i] > 0 ) 
			return (unsigned int)i;
	}
	panic( __FILE__, __LINE__ );
	return INT_MAX; // just to silence compiler
}

/**
  * Notice that all calcuations are performed without actually requiring
  * an array of ranks.
  */
unsigned int MixCovars::rank_sums( 
	double *sums ) {

	const unsigned int N = size();
	const bool PRECALC_SPEARMAN_PRODUCT
		= categoricalIsBinary();
	const unsigned UPPER
		= edge[1].index;
	unsigned int ties = 0;
	unsigned int until = 0;
	float rank = 0;
	double diff;

	// The default pair::operator< automatically uses the first element.

	sort( samples.begin(), samples.end() );

	for(unsigned int i = 0; i < N; i++ ) {

		const unsigned int cat 
			= samples[i].second;

		// The rank of the current sample is i+1 UNLESS we're in
		// the midst of a run of samples with the same value in
		// which case the last computed rank (a mean) is reused.

		if( not ( i < until ) ) {

			until = i + 1;
			while( until < N 
					and ( samples[i].first == samples[until].first ) ) {
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
		diff = rank - meanRank;

		sum_sq_dev += ( diff * diff );

		if( PRECALC_SPEARMAN_PRODUCT ) {
			sum_dev_prod += ( diff * ( edge[UPPER==cat?1:0].meanRank - meanRank ) );
		}

		sums[ cat ] += rank;
	}

	return ties;
}


/**
  * - Establish whether or not Spearman is even sensible (according to
  *   whether or not the categorical variable is binary).
  * - Compute the mean rank used throughout subsequent calcs.
  * - Compute the mean ranks of the categorical values assuming
  *   they -are- binary.
  */
bool MixCovars::complete() {

	if( degenerate() ) return false;

	meanRank = (1.0 + samples.size())/2.0;

	edge[1].index = maxCat(); // ...this is used by kruskal_wallis, but

	assert( edge[1].index < expected_categories );

	// ...the remainder of the edge[] struct is only used if...

	if( categoricalIsBinary() ) {

		// All the following is relevant only to Spearman rho calculation.

		edge[0].index = minCat();
		assert( edge[0].index < expected_categories );

		edge[0].count = category_count[ edge[0].index ];
		edge[0].meanRank = mean_rank_of_ties( 0,    edge[0].count          );

		edge[1].count = category_count[ edge[1].index ];

		if( edge[0].count < arg_min_mixb_count or
			edge[1].count < arg_min_mixb_count )
			return false; // another degeneracy class

		edge[1].meanRank = mean_rank_of_ties( edge[0].count, samples.size());

		assert( (edge[0].count + edge[1].count) == samples.size() );
	}
	return true;
}


/**
 * Calculates the Kruskal-Wallis statistic and optionally a p-value.
 * Importantly, it does it in a one pass iteration over the data.
 */
int MixCovars::kruskal_wallis( struct Statistic *result/*struct CommonStats *cs, struct KruskalWallisStats *ks*/ ) {

	const unsigned int N 
		= size();
	const size_t SIZEOF_SUMS
		= expected_categories * sizeof(double);
	double *rank_sum
		= (double*)alloca( SIZEOF_SUMS );
	memset( rank_sum, 0, SIZEOF_SUMS   );

	double numerator = 0.0;

	result->extra_value[0] = rank_sums( rank_sum );

	for(unsigned int i = 0; i <= edge[1].index; i++ ) {
		if( category_count[i] > 0 ) {
			rank_sum[i] /= category_count[i];
			double delta 
				= rank_sum[i]-meanRank;
			numerator 
				+= ( category_count[i] * delta * delta );
		}
	}

	result->name
		= "Kruskal-Wallis_K";
	result->sample_count
		= N;
	result->value
		= (N-1) * ( numerator / sum_sq_dev );
	result->probability
		= gsl_cdf_chisq_Q( result->value, observed_categories-1 );

	return 0;
}


/**
 * P-value returned is ONE-SIDED. Double it for two-sided.
 */
#ifdef HAVE_MANN_WHITNEY
int MixCovars::mann_whitney( struct CommonStats *cs, struct MannWhitneyStats *ms ) {

	const size_t SIZEOF_SUMS
		= expected_categories * sizeof(double);
	double *rank_sum
		= (double*)alloca( SIZEOF_SUMS );
	memset( rank_sum, 0, SIZEOF_SUMS   );

	ms->ties
		= rank_sums( rank_sum );

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
double MixCovars::spearman_rho() const {

	/**
	  * denominator requires sum of square deviation for
	  * the binary variables which collapses to...
	  */
	double SQD0 = (edge[0].meanRank - meanRank ); SQD0 *= SQD0;
	double SQD1 = (edge[1].meanRank - meanRank ); SQD1 *= SQD1;

	const double SUMSQD
		= edge[0].count * SQD0 + edge[1].count * SQD1;
	/**
	  * And the sum of squared deviation of the numeric covariate's
	  * ranks (sum_sq_dev) was (necessarily) calculated earlier!
	  */
	return sum_dev_prod / sqrt( sum_sq_dev * SUMSQD );
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

#include <iostream>

int main( int argc, char *argv[] ) {

	struct CommonStats common;
	struct KruskalWallisStats kw;
#ifdef HAVE_MANN_WHITNEY
	struct MannWhitneyStats ms;
#endif
	MixCovars accum( MAX_CATEGORY_COUNT );
	char *line = NULL;
	size_t n = 0;
	FILE *fp 
		= argc > 1
		? fopen( argv[1], "r" )
		: stdin;
	while( getline( &line, &n, fp ) > 0 ) {
		unsigned int cat;
		float        num;
		if( line[0] == '#' ) {
			cout << line;
			continue;
		}
		if( 2 == sscanf( line, "%d\t%f\n", &cat, &num ) ) {
			accum.push( num, cat );
		} else {
			cerr << "Failure parsing line: " << line << endl;
		}
	}
	free( line );
	fclose( fp );

	// Always calculate Kruskal-Wallis since it's valid for any groups
	// count >= 2...
	
	if( accum.kruskal_wallis( &common, &kw ) == 0 ) {
		if( accum.categoricalIsBinary() ) {
			printf( "K=%f, p-value=%f, spearman=%f\n", kw.K, common.P, accum.spearman_rho() );
		} else {
			printf( "K=%f, p-value=%f\n", kw.K, common.P );
		}
	} else {
		printf( "error\n" );
	}

#ifdef HAVE_MANN_WHITNEY
	// Do a Mann-Whitney, too, if there are only 2 groups.
	if( accum.categoricalIsBinary() == 2 ) {
		double statistic = cv.mann_whitney_U( &qp );
		printf( "U=%f, p-value=%f (%d ties)\n", statistic, 2.0*qp.prob, qp.kruskal_test.ties );
	}
#endif
	return 0;
}
#endif

