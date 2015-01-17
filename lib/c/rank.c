
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <errno.h>
#include <assert.h>

#include "rank.h"

typedef struct {
	uint32_t off;
	union {
		float    fval;
		uint32_t ival;
	} u;
} pair_t;

typedef const pair_t PAIR_T;

static int _cmp_pairs( const void *pvl, const void *pvr ) {

	PAIR_T *l = (PAIR_T*)pvl;
	PAIR_T *r = (PAIR_T*)pvr;
	if( l->u.fval == r->u.fval ) {
		return 0;
	} else
		return l->u.fval < r->u.fval ? -1 : +1;
}


#ifdef RANK_EMIT_STDERR_WARNINGS
static const char *DEGEN_WARNING 
	= "warning: degenerate case: all %d samples identical\n";
#endif

#if 0
static void _strided_copy( const void *dst, size_t dstride, void *src, size_t sstride, size_t dsize, int n ) {
	const size_t DOFF = dsize*dstride;
	const size_t SOFF = dsize*sstride;
	const uint8_t *d = dst;
	      uint8_t *s = src;
	while( n-- > 0 ) {
		memcpy( d, s, dsize );
		d += DOFF;
		s += SOFF;
	}
}
#endif

void *rank_alloc( int n ) {
	assert( n > 0 );
	return calloc( n, sizeof(pair_t) );
}

void rank_free( void *p ) {
	if( p ) free( p );
}

/**
 * Initializes the pair_t buffer with values from a, possibly discontiguous,
 * memory array. (Discontiguity allowed to accomodate matrices.)
 * The handling of ties is appropriate for use as a prelude to Spearman
 * correlation computation.
 */
int rank_floats( float *base, const unsigned int N, int normalize, void *workspace ) {

	const float NORMALIZER = normalize ? N : 1.0;

	unsigned int i, until = 0;
	pair_t *buf = (pair_t*)workspace;
	float rank, *pf = base;
	int status = 0;

	for(i = 0; i < N; i++ ) {
		buf[i].off = i;
		buf[i].u.fval = *pf++;
	}

	// Sort to determine ranks, and simultaneously check for existance of
	// ties. Absence of ties saves significant work in the sequel.
	
	qsort( buf, N, sizeof(pair_t), _cmp_pairs );

	// We're now done with the -original- floating point values that
	// reside in the pair_t's of buf, so we OVERWRITE those values
	// with their normalized mean rank. These values may then be copied
	// IN ORDER to their final destination.

	for(i = 0; i < N; i++ ) {

		if( ! ( i < until ) ) {

			// Scan forward until we read the first NON-TIE (or
			// the end of data). 
			
			until = i + 1;
			while( until < N 
					&& ( buf[i].u.fval == buf[until].u.fval ) ) {
				status |= RANK_STATUS_TIES;
				until++;
			}

			if( i == 0 && until == N ) {
				status |= RANK_STATUS_CONST;
#ifdef RANK_EMIT_STDERR_WARNINGS
				fprintf( stderr, DEGEN_WARNING, N );
#endif
			}

			// Mean of n integers starting on r:
			// [ (r+0) + (r+1) + (r+2) + ... + (r+(n-1)) ] / n
			// ... = \frac{nr + \sum_{i=0}^{n-1} i}{n} 
			// ... = r + (1/n)[ (n-1)n / 2 ]
			// ... = r + (n-1)/2

			rank = (1.0+i) + ((until-i)-1)/2.0;

			// ...notice if the very next sample -is- different
			// then rank == i+1, the trivial (and usual) case.
		}

		base[ buf[i].off ] = rank / NORMALIZER;
	}

	return status;
}


/**
 * Initializes the pair_t buffer with values from a, possibly discontiguous,
 * memory array. (Discontiguity allowed to accomodate matrices.)
 * The handling of ties is appropriate for use as a prelude to Spearman
 * correlation computation.
 */
int rank_floats_strided( float *base, const unsigned int N, int stride, int normalize, void *workspace ) {

	const float NORMALIZER = normalize ? N : 1.0;

	unsigned int i, until = 0;
	pair_t *buf = (pair_t*)workspace;
	float rank, *pf = base;
	int status = 0;

	assert( stride > 0 );

	for(i = 0; i < N; i++ ) {
		buf[i].off = i;
		buf[i].u.fval = *pf;
		pf += stride;
	}

	// Sort to determine ranks, and simultaneously check for existance of
	// ties. Absence of ties saves significant work in the sequel.
	
	qsort( buf, N, sizeof(pair_t), _cmp_pairs );

	// We're now done with the -original- floating point values that
	// reside in the pair_t's of buf, so we OVERWRITE those values
	// with their normalized mean rank. These values may then be copied
	// IN ORDER to their final destination.

	for(i = 0; i < N; i++ ) {

		if( ! ( i < until ) ) {

			// Scan forward until we read the first NON-TIE (or
			// the end of data). 
			
			until = i + 1;
			while( until < N 
					&& ( buf[i].u.fval == buf[until].u.fval ) ) {
				status |= RANK_STATUS_TIES;
				until++;
			}

			if( i == 0 && until == N ) {
				status |= RANK_STATUS_CONST;
#ifdef RANK_EMIT_STDERR_WARNINGS
				fprintf( stderr, DEGEN_WARNING, N );
#endif
			}

			// Mean of n integers starting on r:
			// [ (r+0) + (r+1) + (r+2) + ... + (r+(n-1)) ] / n
			// ... = \frac{nr + \sum_{i=0}^{n-1} i}{n} 
			// ... = r + (1/n)[ (n-1)n / 2 ]
			// ... = r + (n-1)/2

			rank = (1.0+i) + ((until-i)-1)/2.0;

			// ...notice if the very next sample -is- different
			// then rank == i+1, the trivial (and usual) case.
		}

		base[ stride*(buf[i].off) ] = rank / NORMALIZER;
	}

	return status;
}


#ifdef AUTOUNIT_TEST_RANK

#include <stdio.h>

static int _execR( const char *cmd, int n, float *buf ) {
	float *answer
		= buf + n;
	FILE *fp = popen( cmd, "r" );
	if( fp ) {
		char *line = NULL;
		size_t n = 0;
		int i = 0;
		while( getline( &line, &n, fp ) > 0 ) {
			sscanf( line, "%f\t%f\n", buf+i, answer+i );
			//fprintf( stdout, "%f\t%f\n", buf[i], answer[i] );
			++i;
		}
		if( line )
			free( line );
		pclose( fp );
	} else {
		fprintf( stderr, cmd );
		return -1;
	}
	return 0;
}

static int _test( const char *cmd, int n, float *buf, float *scratch ) {

	float *answer
		= buf + n;

	if( _execR( cmd, n, buf ) )
		return -1;

	rank_floats( buf, n, 0 /* don't normalize */, scratch );

	for(int i = 0; i < n; i++ ) {
		if( buf[i] != answer[i] ) {
			fprintf( stdout, "%f\t != %f\n", buf[i], answer[i] );
			return -1;
		}
	}
	return 0;
}


int main( int argc, char *argv[] ) {

	static char cmd[256];
	const char *CMD_TEMPLATE
		= "R --slave -e 'x<-round(runif(%d,0,%d));"
		"write.table(cbind(x,rank(x)),stdout(),col.names=F,row.names=F,sep=\"\\t\")'";

	const int N
		= argc > 1 ? atoi(argv[1]) : 10;
	const int MAX
		= argc > 2 ? atoi(argv[2]) : 10;
	int ntests
		= argc > 3 ? atoi(argv[3]) : 1;
	float *values
		= calloc( 2*N, sizeof(float) );
	void *workspace = rank_alloc( N );

	sprintf( cmd, CMD_TEMPLATE, N, MAX );

	while( ntests-- > 0 ) {
		if( _test( cmd, N, values, workspace ) )
			break;
	}

	rank_free( workspace );
	free( values );

	return 0;
}

#elif defined(_UNITTEST_RANK_)
// A <- commandArgs( trailingOnly=TRUE );
// N <- if(length(A) > 0 ) as.integer(A[[1]]) else 10;
// f <- if(length(A) > 1 )            A[[2]]  else "foo.tab";
// x <- runif(10)*N
// y <- sample( x, 2*N, replace=TRUE )
// write.table( cbind( y, rank(y) ), f, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
//
// Then use the following bash to verify identical results.
//
// cut -f 1 foo.tab | ./ut_rank $(cat foo.tab | wc -l) | paste foo.tab - 

int main( int argc, char *argv[] ) {

	if( argc > 1 ) {

		const int N = atoi( argv[1] );
		char *line = NULL;
		size_t n = 0;
		int i = 0;
		float *fbuf = (float*)calloc( N, sizeof(float) );
		void *workspace = rank_alloc( N );

		FILE *fp 
			= argc > 2
			? fopen( argv[2], "r" )
			: stdin;

		while( i < N && getline( &line, &n, fp ) > 0 ) {
			if( line[0] == '#' ) {
				fputs( line, stdout );
				continue;
			}
			fbuf[i++] = atof( line );
		}
		free( line );
		fclose( fp );

		rank_floats( fbuf, N, 0 /* don't normalize */, workspace );

		for(i = 0; i < N; i++ ) {
			fprintf( stdout, "%f\n", fbuf[i] );
		}
		rank_free( workspace );
		free( fbuf );
	}
	return EXIT_SUCCESS;
}

#endif

