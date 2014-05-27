
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <alloca.h>

#include <gsl/gsl_cdf.h>

#include "stattest.h"
#include "cat.h"
#include "fisher.h"
#include "min2.h"
#include "bvr.h"

typedef unsigned int count_t;
typedef double       prob_t;

/**
  * This class is stateful in the interest of efficiency
  * because some methods share preliminary computations:
  * 1) statistics require expectation, and 
  * 2) expectation requires marginals
  * (which may also be used incidentally for other purposes), so this 
  * variable captures the progress towards complete calculations
  * to enable implicit just-in-time calculation.
  */
enum State {
	Nothing = 0,
	Marginals,
	Expectation
};


struct CatCovars {

	/**
	 * Buffers are sized to accomodate these dimensions at construction
	 * and are not changeable thereafter.
	 */
	unsigned int ROW_CAPACITY, COL_CAPACITY;

	/**
	  * "decl(ared)_" because the bounds says nothing about the row/col 
	  * content--empty rows/columns may exist! I'm just driving this point
	  * home in these members' names because it has already burned me once!
	  */
	unsigned int decl_rows, decl_cols, decl_cells;

	/**
	  * Count of samples pushed since last cat_clear.
	  */
	count_t sample_count;

	// Calculation state ///////////////////////////////////////////////////

	prob_t minimumExpected;

	/**
	 * This is used strictly by badCellExists
	 */
	count_t REQUIRED_CELL_MINIMUM;

	enum State calculated;

	// Buffers /////////////////////////////////////////////////////////////

	size_t SIZEOF_PROB_BUF;

	prob_t  *expect;

	size_t SIZEOF_COUNT_BUF;

	/**
	  * The counts matrix is stored packed against the head of this
	  * (heap-allocated) buffer. In other words, though this buffer
	  * is -expected- to be larger than required (since the class 
	  * is designed to be iteratively re-used), only the first 
	  * (decl_rows*decl_cols) of the buffer are used for any given 
	  * matrix.
	  */
	count_t *counts;

	/**
	  * Following are derivative data...computed from counts.
	  */
	count_t *rmarg;
	count_t *cmarg;
};


static void _coordsOfOffset(
	struct CatCovars *co, unsigned int off, unsigned int *r, unsigned int *c ) {
	*r = off / co->decl_cols;
	*c = off % co->decl_cols;
}


static count_t _count( struct CatCovars *co, coord_t r, coord_t c ) {
	const int R = co->decl_rows;
	const int C = co->decl_cols;
	assert( r < R && c < C );
	return co->counts[ r*C + c ];
}


static prob_t _expected( struct CatCovars *co, coord_t r, coord_t c ) {
	const int R = co->decl_rows;
	const int C = co->decl_cols;
	assert( r < R && c < C );
	return co->expect[ r*C + c ];
}


static void _recalcSampleCount( struct CatCovars *co ) {
	int n = co->decl_cells;
	co->sample_count = 0;
	while( n-- > 0 ) co->sample_count += co->counts[n];
}


/**
  * This starts searching at the given coordinates and continues linearly
  * through the matrix until it finds a violating cell at which point it
  * turns the linear index into 2D coordinates.
  */
static bool _badCellExists( struct CatCovars *co,
		unsigned int *offset ) {

	const unsigned int N = co->decl_cells;
	unsigned int i = *offset;
	while( i < N ) {
		if( co->counts[i] < co->REQUIRED_CELL_MINIMUM ) {
			*offset = i;
			return true;
		}
		i++;
	}
	return false;
}


/**
  * If after accumulation there are fewer than 2 rows or fewer than 2
  * columns, the table is degenerate. Rows or columns with zero marginal
  * sum are OK! For a contingency table that's bad but not -degenerate-.
  */
static bool _immediatelyDegenerate( struct CatCovars *co ) {
	return ! ( co->decl_rows >= 2 && co->decl_cols >= 2 );
}

#if 0
static bool _eventuallyDegenerate( struct CatCovars *co ) {
	unsigned nonemprows = 0;
	unsigned nonempcols = 0;
	unsigned i;
	for(i = 0; i < co->decl_rows; i++ ) if( co->rmarg[i] > 0 ) nonemprows++;
	for(i = 0; i < co->decl_cols; i++ ) if( co->cmarg[i] > 0 ) nonempcols++;
	return  !  ( nonemprows >= 2 && nonempcols >= 2 );
}
#endif

static void _cullRow( struct CatCovars *co, unsigned int r ) {

	// An over-writing copy suffices since matrix is
	// stored row-major. mergeCols is more intricate.
	
	if( r+1 < co->decl_rows ) { // ...so 0 < rows-U-1
		memmove( 
			co->counts + (r  )*co->decl_cols, 
			co->counts + (r+1)*co->decl_cols, 
			(co->decl_rows - r - 1)*co->decl_cols*sizeof(count_t) );
	} else {
		// ...it's the last row, and there's nothing to move.
	}
	co->decl_rows  -= 1;
	co->decl_cells -= co->decl_cols;

	co->sample_count = 0;
	co->calculated = Nothing; // ...to trigger recomputation.
}


static void _cullCol( struct CatCovars *co, unsigned int c ) {

	const unsigned int N = co->decl_cells;
	// Shrink. Because matrices are stored row-major and calculations
	// throughout this class assume the matrix is packed at the head of
	// the buffer, there is actually a fairly complex shifting of all  
	// but the first c entries that has to occurr...
	// Basically, traverse the matrix linearly (storage order) shifting
	// "left" all entries that are NOT in the c'th column...
	for(unsigned int src = 0, dst = 0 ; src < N; src++ ) {
		if( src % co->decl_cols != c ) {
			co->counts[dst++] = co->counts[src];
		}
	}
	co->decl_cols -= 1;
	co->decl_cells -= co->decl_rows;

	co->sample_count = 0;
	co->calculated = Nothing; // ...to trigger recomputation.
}


/**
  * Following two MUST be called in order before any of
  * the statistics are computed. In other words,
  *    !!!! THESE ARE NOT CALLED IMPLICITLY! !!!!
  * (Statistical methods which -ought- to be const can't be unless
  *  'expect', '[rc]marg' and friends are made mutable->)
  */
static void _calc_marginals( struct CatCovars *co ) {

	const int R = co->decl_rows;
	const int C = co->decl_cols;
	co->sample_count = 0;

	memset( co->rmarg, 0, R*sizeof(count_t) );
	memset( co->cmarg, 0, C*sizeof(count_t) );

	// Calculate the marginals...

	for(unsigned int i = 0; i < R; i++ ) {
		for(unsigned int j = 0; j < C; j++ ) {
			const count_t C = _count( co, i,j);
			co->rmarg[i] += C;
			co->cmarg[j] += C;
			co->sample_count    += C;
		}
	}

	co->calculated = Marginals;
}


/**
 * Calculates the expected probabilities of a contingency table assuming
 * independence of the two variables.
 * "Since the Chi-square and G^2 statistics rest on the normal approximation
 *  to the hypergeometric distribution, these approximations get strained
 *  when there are zero frequency cells or when the expected counts are
 *  'very small'." ...from Yaramakala's Masters Thesis, p.43
 *
 * "...for an independence test to be reliable, all cells have non-zero
 *  expected values and at least 80% of the cells should have expected
 *  values greater than 5." [Everitt 1977, Upton 1978]
 */
static void _calc_expectation( struct CatCovars *co ) {

	const int R = co->decl_rows;
	const int C = co->decl_cols;

	if( co->calculated < Marginals ) 
		_calc_marginals( co );

	if( co->sample_count > 0 ) {	
		co->minimumExpected 
			= (prob_t)(co->rmarg[0] * co->cmarg[0]) 
			/ (prob_t)co->sample_count;
		for(unsigned int i = 0; i < R; i++ ) {
			for(unsigned int j = 0; j < C; j++ ) {
				const prob_t E 
					= (prob_t)(co->rmarg[i] * co->cmarg[j]) 
					/ (prob_t)co->sample_count;
				co->expect[i*co->decl_cols + j] = E;
				if( co->minimumExpected > E ) 
					co->minimumExpected = E;
			}
		}
	}

	co->calculated = Expectation;
}


/**
 * This function cannot fail. Valid results always exist.
 * Returns the indices of the minimal row and columns marginals.
 */
#if 0
static void _minimalMarginals( struct CatCovars *co, unsigned int *rm, unsigned int *cm ) {

	unsigned int MIN_DIM
		= co->decl_rows < co->decl_cols ? co->decl_rows : co->decl_cols;
	unsigned int i, ri = 0, ci = 0;
	count_t mr, mc;

	if( co->calculated < Marginals )
		_calc_marginals( co );

	mr = co->rmarg[ri]; 
	mc = co->cmarg[ci];

	// Iterate over marginals at indices [0,min(rows,cols) )

	for(i = 0; i < MIN_DIM; i++ ) {
		if( mr > co->rmarg[i] ) {
			mr = co->rmarg[i];
			ri = i;
		}
		if( mc > co->cmarg[i] ) {
			mc = co->cmarg[i];
			ci = i;
		}
	}

	// ...continue iteration over [ min(rows,cols), max(rows,cols) ).

	if( co->decl_cols > co->decl_rows )
		for(; i < co->decl_cols; i++ ) {
			if( mc > co->cmarg[i] ) {
				mc = co->cmarg[i];
				ci = i;
			}
		}
	else
	if( co->decl_rows > co->decl_cols )
		for(; i < co->decl_rows; i++ ) {
			if( mr > co->rmarg[i] ) {
				mr = co->rmarg[i];
				ri = i;
			}
		}

	*rm = ri;
	*cm = ci;
}
#endif


/**
 * A pretty printer...primarily for debugging.
 */
#if defined(_UNITTEST_CAT_)
static void dbg_dump( struct CatCovars *co, FILE *fp ) {

	fprintf( fp, "%d(r) x %d(c) (capacity = %dx%d)\n",
		co->decl_rows,
		co->decl_cols,
		co->ROW_CAPACITY,
		co->COL_CAPACITY );

	fprintf( fp, "cell counts...\n" );
	for(unsigned int r = 0; r < co->decl_rows; r++ ) {
		for(unsigned int c = 0; c < co->decl_cols; c++ ) {
			fprintf( fp, "%d\t", _count(co,r,c) );
		}
		fputc( '\n', fp );
	}

	if( co->calculated >= Expectation ) {
		fprintf( fp, "%d total samples, expectation...\n", co->sample_count );
		for(unsigned int r = 0; r < co->decl_rows; r++ ) {
			for(unsigned int c = 0; c < co->decl_cols; c++ ) {
				fprintf( fp, "%.3e\t", _expected( co, r, c) );
			}
			fputc( '\n', fp );
		}
	}
}
#endif


/***************************************************************************
  * Publics
  */

void cat_destroy( void *pv ) {

	if( pv ) {
		struct CatCovars *co = (struct CatCovars *)pv;
		if( co->expect )
			free( co->expect );
		free( pv );
	}
}


/**
  * Pre-allocate a set of working buffers large enough for all anticipated
  * calculations (max feature length) and a struct to wrap them.
  */
void *cat_create( unsigned int rcap, unsigned int ccap ) {

	struct CatCovars *co
		= calloc( 1, sizeof(struct CatCovars) );

	if( co ) {
		co->ROW_CAPACITY = rcap;
		co->COL_CAPACITY = ccap;
		co->SIZEOF_COUNT_BUF
			= (rcap*ccap+rcap+ccap)*sizeof(count_t);
		co->SIZEOF_PROB_BUF
			= rcap*ccap*sizeof(prob_t);

		// Allocate one large buffer and partition it up.
		co->expect
			= calloc( co->SIZEOF_PROB_BUF + co->SIZEOF_COUNT_BUF, sizeof(char) );
		co->counts
			= (count_t*)( co->expect + rcap*ccap );
		co->rmarg
			= co->counts + rcap*ccap;
		co->cmarg
			= co->rmarg + rcap;

		co->REQUIRED_CELL_MINIMUM = 5;

		// If -anything- failed clean up any successes.
		if( (NULL == co->counts) || 
			(NULL == co->expect) ) {
			cat_destroy( co );
			return NULL;
		}
		return co;
	}
	return NULL;
}


void cat_clear( void *pv, coord_t nr, coord_t nc ) {

	struct CatCovars *co = (struct CatCovars *)pv;

	assert( nr <= co->ROW_CAPACITY && nc <= co->COL_CAPACITY );

	// Set dimensions and cell count.

	co->decl_rows  = nr;
	co->decl_cols  = nc;
	co->decl_cells = nr*nc;

	co->sample_count    = 0;
	co->minimumExpected = 0.0;
	co->calculated = Marginals;

	// OPTIMIZATION: Only clearing the part of capacity intended for use!

	memset( co->counts, 0, co->decl_cells*sizeof(count_t) );
	memset( co->rmarg,  0, co->decl_rows*sizeof(count_t) );
	memset( co->cmarg,  0, co->decl_cols*sizeof(count_t) );
	memset( co->expect, 0, co->decl_cells*sizeof(prob_t) );
}


void cat_setMinCellCount( void *pv, unsigned n ) {
	struct CatCovars *co = (struct CatCovars *)pv;
	co->REQUIRED_CELL_MINIMUM = n;
}


/**
 * Only this one is public in order that we may count the samples
 * as they're added...
 */
void cat_push( void *pv, coord_t r, coord_t c ) {

	struct CatCovars *co = (struct CatCovars *)pv;

	assert( r < co->decl_rows );
	assert( c < co->decl_cols );

	co->counts[ r*co->decl_cols + c ] += 1;

	// size and marginals are kept up to date...
	co->sample_count++;
	co->rmarg[r] += 1;
	co->cmarg[c] += 1;
}


size_t cat_size( void *pv ) {
	struct CatCovars *co = (struct CatCovars *)pv;
	return co->sample_count;
}


bool cat_complete( void *pv ) {
	struct CatCovars *co = (struct CatCovars *)pv;
	return ! _immediatelyDegenerate( co );
}


bool cat_is2x2( void *pv ) {
	struct CatCovars *co = (struct CatCovars *)pv;
	return (2 == co->decl_rows) && (2 == co->decl_cols);
}


/**
 * This removes all rows and any columns containing cells with
 * counts LESS THAN OR EQUAL TO REQUIRED_CELL_MINIMUM.
 * 
 * This method assumes we start with a table larger in at least
 * one dimention than 2x2 and halts when/if we reach a 2x2.
 */
unsigned int cat_cullBadCells( void *pv, char *log, int buflen ) {

	struct CatCovars *co = (struct CatCovars *)pv;
	unsigned int culled = 0;
	// BEGIN log maintenance...
	const char * const EOL = log + buflen - 2; // leave room for "+\0".
	const int MIN_LOG_RECORD_LEN = 3; // "R99"
	// ...since tables will never exceed 99 rows or columns.
	char *pc = log;
	// ...END log maintenance.

	unsigned int victim;
	unsigned int linoff = 0; // linear offsets of "bad" cell.

	assert( ! _immediatelyDegenerate( co ) ); // so BOTH dims >= 2

	while( (co->decl_rows > 2 || co->decl_cols > 2) && _badCellExists( co, &linoff ) ) {

		char cullty = '?';
		unsigned int r, c;
		_coordsOfOffset( co, linoff, &r, &c );

		if( co->calculated < Marginals ) 
			_calc_marginals( co );

#ifdef _DEBUG
		if( getenv("SHOW_CULLING") ) {
			fprintf( stderr, 
				"%d at linear offset %d (row %d, col %d) is < %d. R/C marginals are %d/%d\n",
				co->counts[linoff], linoff, 
				r, c, co->REQUIRED_CELL_MINIMUM,
				co->rmarg[r],
				co->cmarg[c] );
		}
#endif

		if( co->rmarg[r] < co->cmarg[c] ) { // Prefer to cull row
			if( co->decl_rows > 2 ) {
				_cullRow( co, (victim = r) );
				cullty = 'R';
			} else 
			if( co->decl_cols > 2 ) {
				_cullCol( co, (victim = c) );
				cullty = 'C';
			} else
				break;
		} else {                    // Prefer to cull column
			if( co->decl_cols > 2 ) {
				_cullCol( co, (victim = c) );
				cullty = 'C';
			} else 
			if( co->decl_rows > 2 ) {
				_cullRow( co, (victim = r) );
				cullty = 'R';
			} else
				break;
		}

		// Log it, if possible...

		if( pc ) {

			pc += sprintf( pc, "%c%d", cullty, victim );

			// Following is a bit of overkill to make sure we can
			// give some indication that we ran out of log space.

			if( EOL - pc < MIN_LOG_RECORD_LEN ) {
				if( EOL - pc > 0 ) strcpy( pc, "+");
				pc = NULL; // no more logging.
				// Really this should never be executed if log buffer
				// is always sufficient.
			}
		}

		// Unlike merging, decimating CAN leave a 0 cell at a lower
		// linear index, so reset is required...
		culled++;
		linoff = 0;
	}

	if( culled > 0 ) 
		_recalcSampleCount( co );

	return culled;
}


#ifdef HAVE_G_TEST
int cat_g( void *pv, struct Statistic *result ) {

	struct CatCovars *co = (struct CatCovars *)pv;
	const int R = co->decl_rows;
	const int C = co->decl_cols;
	double g = 0.0;

	abort(); // unfinished.

	if( co->calculated < Expectation )
		_calc_expectation( co );

	for(unsigned int i = 0; i < R; i++ ) {
		for(unsigned int j = 0; j < C; j++ ) {
			const double O = _count( co, i, j );
			const double E = _expected( co, i, j );
			if( E > 0 && O > 0 ) 
				g += ( O * log( O / E ) );
		}
	}

	g *= 2.0; // TODO: verify!

	result->name
		= "G";
	result->sample_count
		= co->sample_count;
	result->probability
		= gsl_cdf_chisq_Q( g, (R-1)*(C-1) );

	return 0;
}
#endif


int cat_chi_square( void *pv, struct Statistic *result ) {

	struct CatCovars *co = (struct CatCovars *)pv;
	const int R = co->decl_rows;
	const int C = co->decl_cols;

	unsigned int n_empty = 0;
	double chi = 0.0;

	if( co->calculated < Expectation )
		_calc_expectation( co );

	for(unsigned int i = 0; i < R; i++ ) {
		for(unsigned int j = 0; j < C; j++ ) {
			const double O = _count( co, i, j );
			const double E = _expected( co, i, j );
			if( E > 0 && O > 0 ) 
				chi += ((O-E)*(O-E))/E;
			else
			if( O == 0 )
				n_empty++;
		}
	}

	result->name
		= "Chi-square";
	result->sample_count
		= co->sample_count;
	result->probability
		= gsl_cdf_chisq_Q( chi, (R-1)*(C-1) );

	result->extra_value[0] = R;
	result->extra_value[1] = C;
	result->extra_value[2] = co->minimumExpected;
	result->extra_value[3] = n_empty;

	return 0;
}


int cat_fisher_exact( void *pv, struct Statistic *result ) {

	struct CatCovars *co = (struct CatCovars *)pv;

	// OPTIMIZE: does not make sense to call fisher exact without
	// NULL stest_t argument. Keeping this signature for consistency
	// with all other stats, but it might make sense to change the
	// signature of just this function to return P-value (unless I
	// ever added log-likelihood ratio calculation which is what
	// ought to be returned as the function value).

	/**
	  * Given the table:
	  *    a | b
	  *    --+--
	  *    c | d
	  *
	  * ...we need to pass fexact_prob( x, m, n, k )
	  * m == a+c
	  * n == b+d
	  * k == a+b
	  * x == a
	  */
	result->name
		= "Fisher_Exact";
	result->sample_count
		= co->sample_count;
	result->probability
		= fexact_prob( 
			_count(co,0,0), 
			_count(co,0,0) + _count(co,1,0),   // a+c
			_count(co,0,1) + _count(co,1,1),   // b+d
			_count(co,0,0) + _count(co,0,1) ); // a+b

	// TODO: Don't need these now.
	result->extra_value[0] = co->decl_rows;
	result->extra_value[1] = co->decl_cols;

	return 0;
}


/**
  */
double cat_spearman_rho( void *pv ) {

	struct CatCovars *co = (struct CatCovars *)pv;
	double rho = nan("nan");
	const unsigned r0 = _count(co,0,0) + _count(co,0,1);
	const unsigned r1 = _count(co,1,0) + _count(co,1,1);
	const unsigned c0 = _count(co,0,0) + _count(co,1,0);
	const unsigned c1 = _count(co,0,1) + _count(co,1,1);

	// Continuing is pointless if any of these are 0.

	if( r0 > 0 && r1 > 0 && c0 > 0 && c1 > 0 ) {

		const float MU    = (co->sample_count+1.0)/2.0;
		const float r0_mu = (r0+1.0)/2.0;
		const float r1_mu = MEAN_RANK_OF_TIES( r0, co->sample_count );
		const float c0_mu = (c0+1.0)/2.0;
		const float c1_mu = MEAN_RANK_OF_TIES( c0, co->sample_count );

		/**
		  * Following is a heavily compacted form of the Pearson rho,
		  * compacted because thanks to all the ties there are only two
		  * "versions" of (x_i-x_mu) and (y_i-x_mu), and the general
		  * sums in the Pearson quotient just become multiples of these.
		  */
		const double DENOM
			= sqrt(
			  (r0*(r0_mu-MU)*(r0_mu-MU) + r1*(r1_mu-MU)*(r1_mu-MU))
			* (c0*(c0_mu-MU)*(c0_mu-MU) + c1*(c1_mu-MU)*(c1_mu-MU)) );

		if( isnormal(DENOM) ) {
			rho =
				( _count(co,0,0)*(r0_mu-MU)*(c0_mu-MU)
				+ _count(co,0,1)*(r0_mu-MU)*(c1_mu-MU)
				+ _count(co,1,0)*(r1_mu-MU)*(c0_mu-MU)
				+ _count(co,1,1)*(r1_mu-MU)*(c1_mu-MU) )
				/
				DENOM;
		}
		
	}
	return rho;
}


#ifdef _UNITTEST_CAT_

#include <err.h>

int main( int argc, char *argv[] ) {

	if( argc >= 3 ) {

		const int MERGE_LOG_LEN = 128;
		char merge_log[ MERGE_LOG_LEN ]; 

		struct Statistic result;

		const unsigned int MAXROW
			= atoi(argv[1]);
		const unsigned int MAXCOL
			= atoi(argv[2]);

		FILE *fp 
			= argc > 3
			? fopen( argv[3], "r" )
			: stdin;

		int err = 0;
		char *line = NULL;
		size_t n = 0;

		void *accum
			= cat_create( MAXROW, MAXCOL);

		cat_setMinCellCount( accum, 5 );
		cat_clear( accum, MAXROW, MAXCOL );   // ...as per command line.

		while( getline( &line, &n, fp ) > 0 ) {
			int cat[2];
			if( line[0] == '#' ) {
				fputs( line, stdout );
				continue;
			}
			if( 2 == sscanf( line, "%d\t%d\n", cat+0, cat+1 ) ) {
				if( 0 <= cat[0] && 0 <= cat[1] ) {
					cat_push( accum, cat[0], cat[1] );
				}
			} else {
				fprintf( stderr, "Failure parsing line: %s", line );
			}
		}
		free( line );
		fclose( fp );

		puts( "Original..." );
		dbg_dump( accum, stdout );

		if( cat_cullBadCells( accum, merge_log, MERGE_LOG_LEN ) > 0 ) {
			fprintf( stdout, "Culled %s, leaving...\n", merge_log );
			dbg_dump( accum, stdout );
		}

		memset( &result, 0, sizeof(struct Statistic) );
		err = cat_is2x2( accum )
			? cat_fisher_exact( accum, &result )
			: cat_chi_square( accum, &result );

		if( ! err ) {
			printf( "p-value=%f", result.probability );
			if( cat_is2x2( accum ) )
				printf( ", spearman=%f\n", cat_spearman_rho( accum ) );
			else
				fputc( '\n', stdout );
		} else
			printf( "error\n" );

		cat_destroy( accum );

	} else
		err( -1, "%s <max row> <max col> [ <2-column file> ]", argv[0] );
	return 0;
}
#endif

