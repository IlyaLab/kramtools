
#ifdef _UNITTEST_CAT_
#include <iostream>
#else
#include <ostream>
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <alloca.h>

#include <gsl/gsl_cdf.h>

#include "stattest.h"
#include "cat.h"
#include "fisher.h"
#include "min2.h"
#include "bvr.h"

using namespace std;

/**
 * A pretty printer...primarily for debugging.
 */
#if defined(_DEBUG) || defined(_UNITTEST_CAT_)
ostream& operator<<( ostream& os, const CatCovars& ob ) {

	os  << ob.decl_rows << "(r) x " 
		<< ob.decl_cols << "(c) (capacity = " 
		<< ob.row_capacity << "x" << ob.col_capacity
		<< ")" 
		<< endl;

	os << "cell counts..." << endl;
	for(unsigned int r = 0; r < ob.decl_rows; r++ ) {
		for(unsigned int c = 0; c < ob.decl_cols; c++ ) {
			os << ob.count(r, c) << "\t";
		}
		os << endl;
	}

	if( ob.calculated >= CatCovars::Expectation ) {
		os << ob.sample_count
			<< " total samples, expectation..." 
			<< endl;
		for(unsigned int r = 0; r < ob.decl_rows; r++ ) {
			for(unsigned int c = 0; c < ob.decl_cols; c++ ) {
				os << ob.expected(r, c) << "\t";
			}
			os << endl;
		}
	}

	return os;
}
#endif

CatCovars::CatCovars( unsigned int rcap, unsigned int ccap ):
	row_capacity(rcap), 
	col_capacity(ccap), 
	minimumAllowedCellCount(5) {

	counts = new count_t[rcap*ccap];
	rmarg  = new count_t[rcap     ];
	cmarg  = new count_t[     ccap];
	expect = new prob_t [rcap*ccap];

	clear( rcap, ccap );
}

CatCovars::~CatCovars() {
	delete []expect;
	delete []cmarg;
	delete []rmarg;
	delete []counts;
}


void CatCovars::setMinCellCount( unsigned n ) {
	minimumAllowedCellCount = n;
}

void CatCovars::clear( coord_t nr, coord_t nc ) {

	decl_rows = nr;
	decl_cols = nc;

	assert( decl_rows <= row_capacity and decl_cols <= col_capacity );

	// OPTIMIZATION: Only clearing part of capacity intended for use!

	memset( counts, 0, decl_rows*decl_cols*sizeof(count_t) );
	sample_count    = 0;
	minimumExpected = 0;
	memset( rmarg,  0, decl_rows*sizeof(count_t) );
	memset( cmarg,  0, decl_cols*sizeof(count_t) );
	calculated = Marginals;
	// ...as long as they're continuously updated in push!
}


bool CatCovars::complete() {
	return not immediatelyDegenerate();
}


void CatCovars::recalcSampleCount() {
	int n = decl_rows*decl_cols;
	sample_count = 0;
	while( n-- > 0 ) {
		sample_count += counts[n];
	}
}


/**
  * This starts searching at the given coordinates and continues linearly
  * through the matrix until it finds a violating cell at which point it
  * turns the linear index into 2D coordinates.
  */
bool CatCovars::badCellExists( 
		unsigned int *offset ) const {

	const unsigned int N = decl_rows*decl_cols;
	unsigned int i = *offset;
	while( i < N ) {
		if( counts[i] < minimumAllowedCellCount ) {
			*offset = i;
			return true;
		}
		i++;
	}
	return false;
}


/**
 * This function cannot fail. Valid results always exist.
 */
void CatCovars::minimalMarginals( unsigned int *rm, unsigned int *cm ) {

	unsigned int MIN_DIM
		= decl_rows < decl_cols ? decl_rows : decl_cols;
	unsigned int i, ri = 0, ci = 0;
	count_t mr, mc;

	if( calculated < Marginals )
		calc_marginals();

	mr = rmarg[ri]; 
	mc = cmarg[ci];

	// Iterate over marginals at indices [0,min(rows,cols) )

	for(i = 0; i < MIN_DIM; i++ ) {
		if( mr > rmarg[i] ) {
			mr = rmarg[i];
			ri = i;
		}
		if( mc > cmarg[i] ) {
			mc = cmarg[i];
			ci = i;
		}
	}

	// ...continue iteration over [ min(rows,cols), max(rows,cols) ).

	if( decl_cols > decl_rows )
		for(; i < decl_cols; i++ ) {
			if( mc > cmarg[i] ) {
				mc = cmarg[i];
				ci = i;
			}
		}
	else
	if( decl_rows > decl_cols )
		for(; i < decl_rows; i++ ) {
			if( mr > rmarg[i] ) {
				mr = rmarg[i];
				ri = i;
			}
		}

	*rm = ri;
	*cm = ci;
}


/**
 * This removes all rows and any columns containing cells with
 * counts LESS THAN OR EQUAL TO minimumAllowedCellCount.
 * 
 * This method assumes we start with a table larger in at least
 * one dimention than 2x2 and halts when/if we reach a 2x2.
 */
unsigned int CatCovars::cullBadCells( char *log, int buflen ) {

	unsigned int culled = 0;
	// BEGIN log maintenance...
	const char * const EOL = log + buflen - 2; // leave room for "+\0".
	const int MIN_LOG_RECORD_LEN = 3; // "R99"
	// ...since tables will never exceed 99 rows or columns.
	char *pc = log;
	// ...END log maintenance.

	unsigned int victim;
	unsigned int linoff = 0; // linear offsets of "bad" cell.

	assert( not immediatelyDegenerate() ); // so BOTH dims >= 2

	while( (decl_rows > 2 or decl_cols > 2) and badCellExists( &linoff ) ) {

		char cullty = '?';
		unsigned int r, c;
		coordsOfOffset( linoff, &r, &c );

		if( calculated < Marginals ) calc_marginals();

#ifdef _UNITTEST_CAT_
		cout << counts[linoff] << " at linear offset " << linoff 
			<< " (row " << r 
			<< ", col " << c << ") is < " << minimumAllowedCellCount << ". Row/col marginals are "
			<< rmarg[r] << "/" 
			<< cmarg[c] << ", resp." << endl;
#endif

		if( rmarg[r] < cmarg[c] ) { // Prefer to cull row
			if( decl_rows > 2 ) {
				cullRow( (victim = r) );
				cullty = 'R';
			} else 
			if( decl_cols > 2 ) {
				cullCol( (victim = c) );
				cullty = 'C';
			} else
				break;
		} else {                    // Prefer to cull column
			if( decl_cols > 2 ) {
				cullCol( (victim = c) );
				cullty = 'C';
			} else 
			if( decl_rows > 2 ) {
				cullRow( (victim = r) );
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
		recalcSampleCount();

	return culled;
}


bool CatCovars::immediatelyDegenerate() const {
	return not ( decl_rows >= 2 and decl_cols >= 2 );
}


bool CatCovars::eventuallyDegenerate() const {
	unsigned nonemprows = 0;
	unsigned nonempcols = 0;
	for(unsigned i = 0; i < decl_rows; i++ ) if( rmarg[i] > 0 ) nonemprows++;
	for(unsigned i = 0; i < decl_cols; i++ ) if( cmarg[i] > 0 ) nonempcols++;
	return not ( nonemprows >= 2 and nonempcols >= 2 );
}


#ifdef HAVE_MERGE_OPTION
#error This is completely untested and probably broken.
long CatCovars::merge( char *log, int buflen ) {

	unsigned int merges = 0;
	// BEGIN log maintenance...
	const char * const EOL = log + buflen - 2; // leave room for "+\0".
	static const char *ROW_FMT = "R%d,%d";
	static const char *COL_FMT = "C%d,%d";
	const int MIN_LOG_RECORD_LEN = 6; // "R99,99"
	// ...since tables will never exceed 99 rows or columns.
	char *pc = log;
	// ...END log maintenance.

	unsigned int off[2]; // Linear offsets of up to 2 "bad" cells.
	off[0] = 0;          // Begin at the top left.

	while( not is2x2() and badCellExists( off+0 ) ) {

		unsigned int a, b;
		bool MERGE_ROWS = true;

		// Start looking for a 2nd bad cell at the next cell.

		off[1] = off[0]+1;

		// If we find a 2nd, then there are either two rows or two columns
		// that MUST be merged.

		if( badCellExists( off+1 ) ) {

			// ...and we found a 2nd!
			
			unsigned int r[2], c[2];
			coordsOfOffset( off[0], r+0, c+0 );
			coordsOfOffset( off[1], r+1, c+1 );
			const bool DIFFERENT_ROWS = r[0] != r[1];
			const bool DIFFERENT_COLS = c[0] != c[1];

			if( DIFFERENT_ROWS and DIFFERENT_COLS ) {
				// ...so we have a choice of what to merge. Keep matrix
				// as square as possible...
				
				if( decl_rows > decl_cols and decl_rows > 2 ) { 
					//MERGE_ROWS = true; initialized at top
					a = r[0]; b = r[1];
				} else 
				if( decl_cols > 2 ) {
					MERGE_ROWS = false;
					a = c[0]; b = c[1];
				} else {
					a = b = 0; // ...implying failure below.
				}

			} else
			if( DIFFERENT_ROWS ) {
				//MERGE_ROWS = true; initialized at top
				a = r[0];
				b = decl_rows > 2 ? r[1] : a; // b==a implies failure below
			} else {
				assert( DIFFERENT_COLS );
				MERGE_ROWS = false;
				a = c[0]; 
				b = decl_cols > 2 ? c[1] : a; // b==a implies failure below
			}

		} else {

			// ...otherwise, there's only one bad cell. Merge it with a 
			// second row or column, depending on what's more abundant.
			// Then choose the row/column with minimal marginal.

			unsigned int r, c;
			coordsOfOffset( off[0], &r, &c );

			calc_marginals();

			if( decl_rows > decl_cols and decl_rows > 2 ) {
				unsigned int pen, min = min2ui( rmarg, decl_rows, &pen );
				a = r;
				b = (a == min ? pen : min);
				//MERGE_ROWS = true; initialized at top
			} else 
			if( decl_cols > 2 ) {
				unsigned int pen, min = min2ui( cmarg, decl_cols, &pen );
				a = c;
				b = (a == min ? pen : min);
				MERGE_ROWS = false;
			} else {
				a = b = 0; // ...implying failure below.
			}
		}

		if( a == b ) {
			// Note that this might only happen after 1 or more
			// successful merges.
			strcpy( pc, "!" );
			recalcSampleCount();
			return MergeFailure; // ...too few of whatever needs merging.
		}

		if( MERGE_ROWS )
			mergeRows( a, b );
		else
			mergeCols( a, b );

		// Log it, if possible...

		if( pc ) {

			if( MERGE_ROWS )
				pc += sprintf( pc, ROW_FMT, a, b );
			else
				pc += sprintf( pc, COL_FMT, a, b );

			// Following is a bit of overkill to make sure we can
			// give some indication that we ran out of log space.

			if( EOL - pc < MIN_LOG_RECORD_LEN ) {
				if( EOL - pc > 0 ) strcpy( pc, "+");
				pc = NULL; // no more logging.
				// Really this should never be executed if log buffer
				// is always sufficient.
			}
		}

		merges++;

		// TODO(?): Reset off[0]? I don't think so merging vectors
		// cannot create new zero entries and it's always "upper"
		// vectors that are decimated, so should not need to reset.
	}

	if( merges > 0 ) 
		recalcSampleCount();

	// If the log cursor (pc) moved (and we made it to this point in the
	// code), successful merging occurred...
	
	return pc != log ? MergeSuccess : MergeUnnecessary;
}


/**
  * Two operations actually occur:
  * 1) The upper-indexed row's values are added to the lower.
  * 2) The matrix is shrunk row_{i-1} <= row_i for i in [U+1,rows).
  *    Since the matrix is stored row-major, this is just a memmove.
  */
void CatCovars::mergeRows( unsigned int a, unsigned int b ) {

	assert( a != b and a < decl_rows and b < decl_rows );

	const unsigned int L = a < b ? a : b; // Lower index
	const unsigned int U = a < b ? b : a; // Upper index

	for(unsigned int c = 0; c < decl_cols; c++ )
		count(L,c) += count(U,c);

	cullRow( U );
}


void CatCovars::mergeCols( unsigned int a, unsigned int b ) {

	assert( a != b and a < decl_rows and b < decl_rows );

	const unsigned int L = a < b ? a : b; // Lower index
	const unsigned int U = a < b ? b : a; // Upper index

	for(unsigned int r = 0; r < decl_rows; r++ )
		count(r,L) += count(r,U);

	cullCol( U );
}
#endif

void CatCovars::cullRow( unsigned int r ) {

	// An over-writing copy suffices since matrix is
	// stored row-major. mergeCols is more intricate.
	
	if( r+1 < decl_rows ) { // ...so 0 < rows-U-1
		memmove( 
			counts + (r  )*decl_cols, 
			counts + (r+1)*decl_cols, 
			(decl_rows-r-1)*decl_cols*sizeof(count_t) );
	} else {
		// ...it's the last row, and there's nothing to move.
	}
	decl_rows -= 1;

	sample_count = 0;
	calculated = Nothing; // ...to trigger recomputation.
}


void CatCovars::cullCol( unsigned int c ) {

	const unsigned int N = decl_rows*decl_cols;
	// Shrink. Because matrices are stored row-major and calculations
	// throughout this class assume the matrix is packed at the head of
	// the buffer, there is actually a fairly complex shifting of all  
	// but the first c entries that has to occurr...
	// Basically, traverse the matrix linearly (storage order) shifting
	// "left" all entries that are NOT in the c'th column...
	for(unsigned int src = 0, dst = 0 ; src < N; src++ ) {
		if( src % decl_cols != c ) {
			counts[dst++] = counts[src];
		}
	}
	decl_cols -= 1;

	sample_count = 0;
	calculated = Nothing; // ...to trigger recomputation.
}


void CatCovars::calc_marginals() {

	sample_count = 0;

	memset( rmarg, 0, decl_rows*sizeof(count_t) );
	memset( cmarg, 0, decl_cols*sizeof(count_t) );

	// Calculate the marginals...

	for(unsigned int i = 0; i < decl_rows; i++ ) {
		for(unsigned int j = 0; j < decl_cols; j++ ) {
			const count_t C = count(i,j);
			rmarg[i] += C;
			cmarg[j] += C;
			sample_count += C;
		}
	}

	calculated = Marginals;
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
void CatCovars::calc_expectation() {

	if( calculated < Marginals ) 
		calc_marginals();

	if( sample_count > 0 ) {	
		minimumExpected 
			= (prob_t)(rmarg[0] * cmarg[0]) 
			/ (prob_t)sample_count;
		for(unsigned int i = 0; i < decl_rows; i++ ) {
			for(unsigned int j = 0; j < decl_cols; j++ ) {
				const prob_t E 
					= (prob_t)(rmarg[i] * cmarg[j]) 
					/ (prob_t)sample_count;
				expect[i*decl_cols + j] = E;
				if( minimumExpected > E ) 
					minimumExpected = E;
			}
		}
	}

	calculated = Expectation;
}


#ifdef HAVE_G_TEST
int CatCovars::g( struct CovarsSummary *cs ) {

	double g = 0.0;

	abort(); // unfinished.
	if( calculated < Expectation )
		calc_expectation();

	for(unsigned int i = 0; i < decl_rows; i++ ) {
		for(unsigned int j = 0; j < decl_cols; j++ ) {
			const double O = count( i, j );
			const double E = expected( i, j );
			if( E > 0 and O > 0 ) 
				g += ( O * log( O / E ) );
		}
	}

	g *= 2.0; // TODO: verify!

	if( qp ) {
		qp->prob = gsl_cdf_chisq_Q( g, (decl_rows-1)*(decl_cols-1) );
	}

	return 0;
}
#endif

int CatCovars::chi_square( struct CommonStats *cs, struct ChiSquareStats *qs ) {

	unsigned int n_empty = 0;
	double chi = 0.0;

	if( calculated < Expectation )
		calc_expectation();

	for(unsigned int i = 0; i < decl_rows; i++ ) {
		for(unsigned int j = 0; j < decl_cols; j++ ) {
			const double O = count( i, j );
			const double E = expected( i, j );
			if( E > 0 and O > 0 ) 
				chi += ((O-E)*(O-E))/E;
			else
			if( O == 0 )
				n_empty++;
		}
	}

	cs->P = gsl_cdf_chisq_Q( chi, (decl_rows-1)*(decl_cols-1) );

	qs->rows = decl_rows;
	qs->cols = decl_cols;
	qs->min_expected_cell_count = minimumExpected;
	qs->empty_cells             = n_empty;

	return 0;
}


int CatCovars::fisher_exact( struct CommonStats *cs, struct FisherExactStats *fs ) const {

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
	cs->P = fexact_prob( 
			count(0,0), 
			count(0,0) + count(1,0),   // a+c
			count(0,1) + count(1,1),   // b+d
			count(0,0) + count(0,1) ); // a+b
	fs->rows = decl_rows;
	fs->cols = decl_cols;

	return 0;
}


/**
  */
double CatCovars::spearman_rho() const {

	double rho = nan("nan");
	const unsigned r0 = count(0,0) + count(0,1);
	const unsigned r1 = count(1,0) + count(1,1);
	const unsigned c0 = count(0,0) + count(1,0);
	const unsigned c1 = count(0,1) + count(1,1);

	// Continuing is pointless if any of these are 0.

	if( r0 > 0 and r1 > 0 and c0 > 0 and c1 > 0 ) {

		const float MU    = (sample_count+1.0)/2.0;
		const float r0_mu = (r0+1.0)/2.0;
		const float r1_mu = mean_rank_of_ties( r0, sample_count );
		const float c0_mu = (c0+1.0)/2.0;
		const float c1_mu = mean_rank_of_ties( c0, sample_count );

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
				( count(0,0)*(r0_mu-MU)*(c0_mu-MU)
				+ count(0,1)*(r0_mu-MU)*(c1_mu-MU)
				+ count(1,0)*(r1_mu-MU)*(c0_mu-MU)
				+ count(1,1)*(r1_mu-MU)*(c1_mu-MU) )
				/
				DENOM;
		}
		
	}
	return rho;
}


#ifdef _UNITTEST_CAT_

/**
 */

#include <cstdio>
#include <cstdlib>

int main( int argc, char *argv[] ) {

	if( argc > 2 ) {

		struct CommonStats common;
		struct ChiSquareStats chisq;
		struct FisherExactStats fisherx;

		const int MERGE_LOG_LEN = 128;
		char merge_log[ MERGE_LOG_LEN ]; 

		const unsigned int MAXROW
			= atoi(argv[1]);
		const unsigned int MAXCOL
			= atoi(argv[2]);

		char *line = NULL;
		size_t n = 0;
		FILE *fp 
			= argc > 3
			? fopen( argv[3], "r" )
			: stdin;
		int err = 0;

		CatCovars table(2*MAXROW,2*MAXCOL); // ...arbitrarily.
		table.setMinCellCount( 5 );
		table.clear( MAXROW, MAXCOL );   // ...as per command line.

		while( getline( &line, &n, fp ) > 0 ) {
			int cat[2];
			if( line[0] == '#' ) {
				cout << line;
				continue;
			}
			if( 2 == sscanf( line, "%d\t%d\n", cat+0, cat+1 ) ) {
				if( 0 <= cat[0] and 0 <= cat[1] ) {
					table.push( cat[0], cat[1] );
				}
			} else {
				cerr << "Failure parsing line: " << line << endl;
			}
		}
		free( line );
		fclose( fp );

		cout << "Original..." << endl;
		cout << table;

		if( table.cullBadCells( merge_log, MERGE_LOG_LEN ) > 0 ) {
			cout << "Culled " << merge_log << ", leaving..." << endl;
			cout << table;
		}

		table.calc_marginals();
		table.calc_expectation();

		err = table.is2x2()
			? table.fisher_exact( &common, &fisherx )
			: table.chi_square( &common, &chisq );

		if( ! err ) {
			if( table.is2x2() )
				printf( "p-value=%f, spearman=%f\n", common.P, table.spearman_rho() );
			else
				printf( "p-value=%f\n", common.P );
		} else
			printf( "error\n" );
	} else
		fprintf( stderr, "%s <max row> <max col> [ <2-column file> ]\n", argv[0] );
	return 0;
}
#endif

