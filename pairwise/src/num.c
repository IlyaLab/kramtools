
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <alloca.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_float.h>

#include "rank.h"
#include "stattest.h"
#include "num.h"

typedef float con_t;

struct ConCovars {

	/**
	  * Maximum number of samples
	  */
	int SAMPLE_CAPACITY;

	/**
	  * Count of samples pushed since last con_clear.
	  */
	int sample_count;

	/**
	  * The total amount of malloc'ed space (for the purposes
	  * of fast clearing with memset).
	  */
	size_t SIZEOF_BUFFERS;

	/**
	  * Two buffers
	  */
	con_t *l, *r;

	void *rank_scratch;
};

#if defined(_UNITTEST_NUM_)
static void dbg_dump( struct ConCovars *co, FILE *fp ) {
	for(unsigned int i = 0; i < co->sample_count; i++ )
		fprintf( fp, "%f\t%f\n", co->l[i], co->r[i] );
}
#endif


void con_destroy( void *pv ) {

	if( pv ) {
		struct ConCovars *co = (struct ConCovars *)pv;
		if( co->rank_scratch )
			rank_free( co->rank_scratch );
		if( co->l )
			free( co->l );
		free( pv );
	}
}


/**
  * Pre-allocate a set of working buffers large enough for all anticipated
  * calculations (max feature length) and a struct to wrap them.
  */
void *con_create( unsigned int cap ) {
	struct ConCovars *co
		= calloc( 1, sizeof(struct ConCovars) );
	if( co ) {
		co->SAMPLE_CAPACITY = cap;
		co->SIZEOF_BUFFERS  = 2*cap*sizeof(con_t);
		// Allocate one large buffer and partition it up.
		co->l = calloc( co->SIZEOF_BUFFERS, sizeof(char) );
		co->r = co->l + cap;
		co->rank_scratch = rank_alloc( cap );
		// If -anything- failed clean up any successes.
		if( (NULL == co->l) || 
			(NULL == co->rank_scratch) ) {
			con_destroy( co );
			return NULL;
		}
		return co;
	}
	return NULL;
}


void con_clear( void *pv ) {
	struct ConCovars *co = (struct ConCovars *)pv;
	co->sample_count = 0;
	memset( co->l, 0, co->SIZEOF_BUFFERS );
}


void con_push( void *pv, float n1, float n2 ) {
	struct ConCovars *co = (struct ConCovars *)pv;
	const int i = co->sample_count++;
	assert( i <= co->SAMPLE_CAPACITY );
	co->l[i] = n1;
	co->r[i] = n2;
}


size_t con_size( void *pv ) {
	return ((struct ConCovars *)pv)->sample_count;
}


bool con_complete( void *pv ) {
	return ((struct ConCovars *)pv)->sample_count > 2;
	// ...otherwise p-value computation will fail.
}


/**
 */
int con_spearman_correlation( void *pv, struct Statistic *result ) {

	struct ConCovars *co = (struct ConCovars *)pv;
	const int N = co->sample_count;

	assert( N > 2 );
	assert( NULL != co->rank_scratch );

	const int rinfo1 
		= rank_floats( co->l, N, 0, co->rank_scratch );
	const int rinfo2 
		= rank_floats( co->r, N, 0, co->rank_scratch );

	if( RANK_STATUS_CONST & rinfo1 ) // vectors were in fact constant!
		result->extra_value[0] = N-1;
	if( RANK_STATUS_CONST & rinfo2 )
		result->extra_value[1] = N-1;

	{
		const double rho 
			= gsl_stats_float_correlation( co->l, 1, co->r, 1, N );

		/**
		 * P-value computation for the correlation.
		 */

#ifdef HAVE_FISHER_TRANSFORM
		const double FisherTransform 
			//= 0.5 * log( (1.+rho) / (1.-rho) );
			= atanh(rho);
		// ...absolute value to simplify CDF use below, since the
		// transformation is symmetric.
		const double z
			= sqrt( (N - 3.0) / 1.06 ) * FisherTransform;
		// ...z ~ N(0,1) under null hyp of statistical independence.
		result->name
			= "Spearman_rho,Fisher_transform";
		result->probability = gsl_cdf_ugaussian_Q( fabs(z) );
#else
		const double t 
			=  fabs( rho*sqrt((N-2.0)/(1.0-rho*rho)) );
		// ...abs so that I can always test the upper tail.
		// x2 below to make it a two-tailed test. (t-distribution
		// is symmetric).
		result->name
			= "Spearman_rho,t-distribution";
		result->probability = 2* gsl_cdf_tdist_Q(t,N-2.0);
#endif
		result->value = rho;
	}
	result->sample_count = N;

	return 0;
}


#ifdef HAVE_SCALAR_PEARSON
int con_pearson_correlation( void *pv, struct Statistic *result ) {
	struct ConCovars *co = (struct ConCovars *)pv;
	const int N = co->sample_count;
	ps->rho = gsl_stats_float_correlation( co->l, 1, co->r, 1, N );
	return 0;
}
#endif


#ifdef _UNITTEST_NUM_

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

	if( argc >= 2 ) {

		struct Statistic result;
		struct ConCovars *accum
			= con_create( atoi(argv[1]) );
		FILE *fp 
			= argc > 2
			? fopen( argv[2], "r" )
			: stdin;
		char *line = NULL;
		size_t n = 0;

		con_clear( accum );

		while( getline( &line, &n, fp ) > 0 ) {
			float l, r;
			if( line[0] == '#' ) {
				fputs( line, stdout );
				continue;
			}
			if( 2 == sscanf( line, "%f\t%f\n", &l, &r ) ) {
				con_push( accum, l, r );
			} else {
				fprintf( stderr, "Failure parsing line: %s", line );
			}
		}
		free( line );
		fclose( fp );

		memset( &result, 0, sizeof(struct Statistic) );
#ifdef HAVE_SCALAR_PEARSON
		con_pearson_correlation( accum, &result );
		printf( "pearson=%f\n", result.value );
#endif
		memset( &result, 0, sizeof(struct Statistic) );
		if( con_spearman_correlation( accum, &result ) == 0 ) {
			printf( "p-value=%f, spearman=%f\n", result.probability, result.value );
		} else
			printf( "error\n" );

		con_destroy( accum );

	} else
		err( -1, "%s <sample count> [ <input file> ]", argv[0] );

	return 0;
}
#endif

