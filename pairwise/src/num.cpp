
#include <ostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <alloca.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_float.h>

#include "rank.h"
#include "stattest.h"
#include "num.h"
#include "fp.h"

using namespace std;

#if defined(_DEBUG) || defined(_UNITTEST_NUM_)
ostream& operator<<( ostream& os, const NumCovars& ob ) {

	for(unsigned int i = 0; i < ob.size(); i++ ) {
		os << ob.l[i] << "\t" << ob.r[i] << endl;
	}
	return os;
}
#endif


NumCovars::NumCovars() :
	rank_scratch(NULL) {
}


NumCovars::~NumCovars() {
	if( rank_scratch )
		rank_free( rank_scratch );
}


void NumCovars::reserve( unsigned int n ) {
	l.reserve( n );
	r.reserve( n );
	rank_scratch = rank_alloc( l.capacity() );
}


void NumCovars::clear() {
	l.clear();
	r.clear();
}


bool NumCovars::complete() {
	return size() > 2; // since p-value computation will fail otherwise
}

/**
 */
int NumCovars::spearman_correlation( struct Statistic *result
		/* struct CommonStats *cs, struct SpearmanStats *ss */ ) {

	const int N
		= size();

	assert( N > 2 );
	assert( NULL != rank_scratch );

	const int rinfo1 
		= rank_floats( l.data(), N, 0, rank_scratch );
	const int rinfo2 
		= rank_floats( r.data(), N, 0, rank_scratch );

	if( RANK_STATUS_CONST & rinfo1 ) // vectors were in fact constant!
		result->extra[0] = N-1;
	if( RANK_STATUS_CONST & rinfo2 )
		result->extra[1] = N-1;

	{
		const double rho 
			= gsl_stats_float_correlation( l.data(), 1, r.data(), 1, N );

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
			= "Spearman rho, Fisher transform";
		result->probability = gsl_cdf_ugaussian_Q( fabs(z) );
#else
		const double t 
			=  fabs( rho*sqrt((N-2.0)/(1.0-rho*rho)) );
		// ...abs so that I can always test the upper tail.
		// x2 below to make it a two-tailed test. (t-distribution
		// is symmetric).
		result->name
			= "Spearman rho, t-dist approx.";
		result->probability = 2* gsl_cdf_tdist_Q(t,N-2.0);
#endif
		result->value = rho;
	}

	return 0;
}


#ifdef HAVE_SCALAR_PEARSON
int NumCovars::pearson_correlation( struct CommonStats *cs, struct PearsonStats *ps ) {
	ps->rho = gsl_stats_float_correlation( l.data(), 1, r.data(), 1, size() );
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

#include <iostream>

int main( int argc, char *argv[] ) {

	struct CommonStats common;
	struct SpearmanStats sp;
	NumCovars cv;
	double statistic;
	char *line = NULL;
	size_t n = 0;
	FILE *fp 
		= argc > 1
		? fopen( argv[1], "r" )
		: stdin;
	while( getline( &line, &n, fp ) > 0 ) {
		float l, r;
		if( line[0] == '#' ) {
			cout << line;
			continue;
		}
		if( 2 == sscanf( line, "%f\t%f\n", &l, &r ) ) {
			cv.push( l, r );
		} else {
			cerr << "Failure parsing line: " << line << endl;
		}
	}
	free( line );
	fclose( fp );

#ifdef HAVE_SCALAR_PEARSON
	statistic = cv.pearson_correlation( &qp );
	printf( "pearson=%f\n", statistic );
#endif
	if( cv.spearman_correlation( &common, &sp ) == 0 ) {
		printf( "p-value=%f, spearman=%f\n", common.P, sp.rho );
	} else
		printf( "error\n" );

	return 0;
}
#endif

