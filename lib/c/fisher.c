
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fisher.h"

#include "gsl/gsl_randist.h"

static double *_buf = NULL;
static unsigned int _allocn = 0;

/**
  * Can't just use _allocn as an indicator since explicit reserve/releases
  * can cycle _allocn between 0 and non-zero, but once atexit is called
  * the registration is permanent.
  */
static int _callback_registered = 0;

static void _fexact_free( void ) {
	if( _buf ) free( _buf );
}

int fexact_reserve( unsigned int n ) {
	_buf = (double*)realloc( _buf, n * sizeof(double) );
	if( _buf ) {
		if( ! _callback_registered ) {
			_callback_registered = 1;
			atexit( _fexact_free );
		}
		_allocn = n;
		return 0;
	} else
		return -1;
}


void fexact_release() {
	if( _buf ) free( _buf );
	_buf = NULL;
	_allocn = 0;
}


double fexact_prob( unsigned int x, unsigned int m, unsigned int n, unsigned int k ) {

	// support is [lo,hi]...
	const unsigned int lo = 0 < k-n ? k-n : 0;
	const unsigned int hi = k < m   ? k   : m;
	const unsigned int N = hi-lo+1;

	double lim, sum = 0.0, max = 0.0;
	int i;

	if( _allocn < N ) {
		_buf = (double*)realloc( _buf, N * sizeof(double) );
		if( NULL == _buf ) {
			fprintf( stderr, "error: failed realloc( %p, %ld ) @ %s:%d\n",
				_buf, N*sizeof(double), __FILE__, __LINE__ );
			return nan("NaN");
		}
		// If this is first allocation, register for cleanup.
		if( ! _callback_registered ) {
			_callback_registered = 1;
			atexit( _fexact_free );
		}
		_allocn = N;
	}

	// Begin calculation.

	i = N;
		_buf[0] 
			= max 
			= log( gsl_ran_hypergeometric_pdf( lo+0, m, n, k ) );
	while( i-- > 1 ) {
		const double d 
			= log( gsl_ran_hypergeometric_pdf( lo+i, m, n, k ) );
		if( max < d )
			max = d;
		_buf[i] = d;
	}
	i = N; while( i-- > 0 ) {
		_buf[i] = exp( _buf[i] - max );
		sum += _buf[i];
	}	
	i = N; while( i-- > 0 ) {
		_buf[i] /= sum;
	}

	lim = _buf[x-lo]*(1.0+1e-7);

	// Sum up the probabilities

	sum = 0.0;
	i = N; while( i-- > 0 ) {
		if( _buf[i] <= lim ) sum += _buf[i];
	}

	return sum;
}


#ifdef _UNITTEST_FISHER_
/**
  * p(k) =  C(n_1, k) C(n_2, t - k) / C(n_1 + n_2, t)
  *
  * ...where C(a,b) = a!/(b!(a-b)!) and t <= n_1 + n_2. 
  * The domain of k is max(0,t-n_2), ..., min(t,n_1).
  *
  *   gsl_ran_hypergeometric_pdf( k,  n1,  n2,  t)
  * This function computes the probability p(k) of obtaining k from a 
  * hypergeometric distribution with parameters n1, n2, t, using the 
  * formula given above.
  *
  * argv[1] | argv[2]      a | b
  * --------+--------  ==  --+--
  * argv[3] | argv[4]      c | d
  */
int main( int argc, char *argv[] ) {

	if( argc > 0 ) {

		char *line = NULL;
		size_t n = 0;
		FILE *fp 
			= argc > 1
			? fopen( argv[1], "r" )
			: stdin;
		while( getline( &line, &n, fp ) > 0 ) {
			unsigned int a, b, c, d;
			if( line[0] == '#' ) {
				fputs( line, stdout );
				continue;
			}
			if( 4 == sscanf( line, "%d\t%d\t%d\t%d", &a, &b, &c, &d ) ) {
				printf( "%f\n",
					fexact_prob( a, a+c /*m*/, b+d /*n*/, a+b /*k*/ )  );
			}
		}
		free( line );
		fclose( fp );

	} else {
		fprintf( stderr, "%s [ <2-column file> ]\n", argv[0] );
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
#endif

