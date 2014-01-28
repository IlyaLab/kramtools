
/**
  * This code supports an extremely flexible (too flexible?) method of
  * generating a sequence of pairs from a range [0,<limit>). It was
  * motivated by need to iterate over subsets of pairs of rows from a
  * matrix, but nothing in here is specific to matrices.
  * It is generally used as:
  *
  * iter_parse( <outer spec>, <inner spec> )
  * if( iter_begin( <limit>, &outer_index, &inner_index ) ) {
  * 	do {
  *			// use outer_index, inner_index
  *		} while( iter_next( &o, &i ) );
  * }
  *
  * The "inner" and "outer" range specifications are extremely flexible.
  * It supports the canonical sort of nested loops like:
  *
  * for(i = 0; i < N-1; i++ ) {
  *   for(j = i+1; j < N; j++ ) {
  *     // do something with i,j
  *   }
  * }
  * ...but
  * 1) the increments need not be one
  * 2) the "inner" loop can, but need not, be relative to the outer
  * 3) and SEQUENCES of ranges are supported
  *    instead of [0,N-1) for the outer loop it could be
  *    [1,3),[10,15),[19,N-1) and the inner loop could be relative
  *    to each of those ranges or an absolute range
  * In other words -almost- arbitrary cross-products are supported.
  * The only absolute requirement is that ranges are internally ascending,
  * though adjacent ranges need not be. e.g. "10:11,9:10" is equivalent
  * to {10,9}
  * 
  * That is, after calling the iter_next method
  * returns PAIRS in [0,<limit>), subject to a specification given as a
  * string to iter_parse.
  *
  * As it is currently implemented (specifically, its reliance on statics
  * and the absence of a reinitialization method) it is intended to be used
  * ONCE in the context of a process. That is command line arguments are 
  * parsed ONCE, the specified iteration executed ONCE, then the process 
  * exits. 
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <regex.h>
#include <assert.h>
#include "iter.h"

#define _minimum(x,y) ((x)<(y)?(x):(y))

typedef struct _interval {
	int start;// ...must always be a legitimate offset in [0,_N).
	int stop; // INVALID only allowed in stop fields.
#define INVALID (-1)
#define is_valid(x) (0<=x)
} interval_t;

typedef struct _range {

	interval_t bounds;

	/**
	  * signed delta's are theoretically possible, but not currently 
	  * supported. iter_next is all geared to increasing indices.
	  * The regexp definition precludes negative numbers.
	  */
	unsigned delta;

	/**
	  * The isrel flag always applies to bounds.start; it applies to bounds.stop
	  * if and only if is_valid(stop). An end-of-data indicator flag
	  * in a relative range implies the -absolute- end of data. What else
	  * could end-of-data mean?
	  */
	bool isrel;

} range_t;

typedef const range_t RANGE_T;

static void _init_range( range_t *r ) {
	memset( r, 0, sizeof(range_t) );
	r->bounds.stop = INVALID; // end of data, whatever that is.
	r->delta       = 1;
}

typedef struct {
	range_t *range;
	int      count;
} range_array_t;

/**
  * private variables
  */
static int _initialized = 0;
static regex_t range_expr;
static int _N = 0; // set by iter_begin.
static range_array_t _outer;
static range_array_t _inner;

#ifdef _DEBUG
static void _iter_dump_range( const char *prefix, 
		const char *which,
		int i, const range_t *r, FILE *fp ) {

	fprintf( fp, "%s%s %02d: %s [%d,%d) by %d\n",
		prefix,
		which,
		i,
		r->isrel ? "REL" : "ABS",
		r->bounds.start, r->bounds.stop,
		r->delta );
}


void iter_dump( const char *prefix, FILE *fp ) {

	int i;
	fprintf( fp, "%s%d outer and %d inner ranges\n",
		prefix,
		_outer.count,
		_inner.count );

	for(i = 0; i < _outer.count; i++ ) {
		_iter_dump_range( prefix, "O", i, _outer.range + i, fp );
	}
	for(i = 0; i < _inner.count; i++ ) {
		_iter_dump_range( prefix, "I", i, _inner.range + i, fp );
	}
}
#endif


static void _iter_fini() {
	regfree( &range_expr );
	if( _inner.range ) free( _inner.range );
	if( _outer.range ) free( _outer.range );
}


/**
  * Current setup for implicitl initialization and finalization.
  * It's done just-in-time in iter_parse. May want to change this.
  */
static int _iter_init() {

	memset( &_inner, 0, sizeof(_inner) );
	memset( &_outer, 0, sizeof(_outer) );

	int regerr = regcomp( 
			&range_expr, 
			"^(\\+[0-9]+|[0-9]*)(:[0-9]*){0,2}(,|$)", 
			REG_EXTENDED | REG_NEWLINE );
	if( regerr ) {
		const size_t required 
			= regerror( regerr, &range_expr, NULL, 0 );
		char *buf 
			= (char *)alloca( required );
		regerror( regerr, &range_expr, buf, required );
		fprintf( stderr, "%s\n", buf );
		return -__LINE__;
	}

	atexit( _iter_fini );
	_initialized = 1;
	return 0;
}


static int _count( const char *s, char c ) {
	int n = 0;
	while( *s ) if( *s++ == c ) n++;
	return n;
}


/**
  */
static int _parse_one( const char *argv, bool inner ) {

	range_array_t * const loop
		= inner
		? &_inner
		: &_outer;
	int i;
	regmatch_t match[1];
	const char *pc;

	if( NULL == argv ) {
		return -__LINE__;
	}

	/**
	 * Interpret the field specification.
	 * First the number of discontiguous ranges...
	 */
	loop->count = _count( argv, ',' ) + 1;

	loop->range 
		= (range_t*)calloc( 
			loop->count+1 /* + a guard */, 
			sizeof(range_t) ); 

	/**
	 * ...then actually parse the ranges.
	 */
	pc = argv;
	for(i = 0; i < loop->count; i++ ) {

		char *sep = NULL;
		range_t *r = loop->range + i;
		_init_range( r );

		if( regexec( &range_expr, pc, 1, match, 0 ) || match[0].rm_so < 0 ) 

			break; // ...we're done

		else { // We have another range expression...

			char rng[16];
			const int L = match[0].rm_eo - match[0].rm_so;

			// Snip the range out of the command line into a working buffer.
			memcpy( rng, pc, L ); rng[L] = 0;

			if( '+' == *pc ) {
				if( inner )
					r->isrel = true;
				else
					fprintf( stderr, 
						"warning: ignoring the '+' on '%s'."
						"(Outer loop cannot be relative.)\n", rng );
			}

			switch( _count( rng, ':' ) ) {

			// Careful below to never initialize using empty fields;
			// LEAVE ANY UNSPECIFIED FIELD IN THE DEFAULT STATE  
			// ESTABLISHED IN _init_range.

			// Following generally works backwards through the range
			// expression, snipping of suffixes...

			case 2:
				sep = strrchr( rng, ':' );
				if( strlen(sep+1) > 0 ) {
					r->delta = (unsigned)atoi( sep+1 );
					// The regexp precludes specification of negative 
					// numbers, so the cast is safe.
				}
				*sep = '\0';
				// fallthrough
			case 1:
				sep = strrchr( rng, ':' );
				if( strlen(sep+1) > 0 ) {
					r->bounds.stop = atoi( sep+1 );
				}
				*sep = '\0';
				// fallthrough
			case 0:
				sep // it's not actually a separator here, just reusing it
					= r->isrel ? rng+1 : rng;
				if( strlen( sep ) > 0 )
					r->bounds.start = atoi( sep );
				break;

			default:
				fprintf( stderr, "error: bad format: '%s'\n", rng );
				return -__LINE__;
			}

#ifdef _UNITTEST_ITER_
			fprintf( stdout, "_parse_one: %s: %02d: %s [%d, %d), delta=%d\n", 
				inner ? "inner":"outer", 
				i, 
				r->isrel? "REL" : "ABS",
				r->bounds.start, r->bounds.stop, 
				r->delta );
#endif
			// Verify that the range bounds, if they're explicit, are sensible.

			static const char *ERRMSG
				= "error: stop (%d) <= start (%d) (%s)\n";
			if( is_valid(r->bounds.stop) && r->bounds.stop <= r->bounds.start ) {
				fprintf( stderr, ERRMSG, r->bounds.stop, r->bounds.start, "absolute" );
				return -__LINE__;
			}
			pc += L;
		}
	}
	return i < loop->count ? -__LINE__ : 0 /* success */;
}


/**
  * The following 4 variables along with the two that the client code
  * maintains FULLY SPECIFY the state of iteration.
  */
static int rt_outer_rng; // index into _outer.range
static int rt_outer_end;

static int rt_inner_rng; // index into _inner.range
static int rt_inner_end;

/**
  * Invariants:
  * rt_???er_end should NEVER exceed _N
  */

/**
  * Update the rt_inner_end and the starting value *offi.
  * Assumes that rt_inner_rng has already been set up.
  * Returns truth of, "the inner iteration range is non-empty."
  */
static bool _begin_inner( const int OUTER, int *offi ) {

	RANGE_T *pI
		= _inner.range + rt_inner_rng;
	const int START
		= pI->isrel ? OUTER : 0;

	assert( rt_inner_rng < _inner.count );

	rt_inner_end
		= is_valid( pI->bounds.stop )
		? _minimum( START + pI->bounds.stop, _N )
		: _N;

	assert( 0 <= pI->bounds.start );

	*offi = START + pI->bounds.start;

	return *offi < rt_inner_end;
}


/***************************************************************************
  * Public API
  */

int iter_parse( const char *outer, const char *inner ) {

	if( ! _initialized ) {
		if( _iter_init() ) return -3;
	}

	if( _parse_one( outer, false ) ) {
		return -1;
	}
	if( _parse_one( inner, true  ) ) {
		return -2;
	}
	return 0;
}


bool iter_begin( int limit, int *offo, int *offi ) {

	RANGE_T *pO
		= &( _outer.range[0] );

	_N = limit;

	// Range indices ALWAYS start at 0.

	rt_outer_rng = 0;
	rt_inner_rng = 0;

	// Outer end is trivial to determine...
	// Outer is NEVER relative.

	rt_outer_end 
		= is_valid( pO->bounds.stop )
		? _minimum( pO->bounds.stop, _N )
		: _N;

	*offo = _outer.range[ 0 ].bounds.start;

	// ...but inner end is less so.

	if( _begin_inner( *offo, offi ) ) {
		return *offo == *offi ? iter_next( offo, offi ) : true;
	}
	return false;
}


/**
  * Inner and outer loop controls always iterate through their
  * ranges exhaustively: outer once, inner as many times as outer
  * dictates.
  * This loop must NEVER return with *offo == *offi.
  */
bool iter_next( int *offo, int *offi ) {

	bool inner_range_changed, outer_range_changed;

	do {

_entry:
		inner_range_changed = false;
		outer_range_changed = false;

		// There are four levels of incrementing to do:
		// 1) inner offset
		// 2) inner range
		// 3) outer offset
		// 4) outer range
		// ...in that order, so first we determine all operations 
		// that will be entailed by incrementing the inner offset.

		*offi += _inner.range[ rt_inner_rng ].delta;

		if( ! ( *offi < rt_inner_end ) ) {
			inner_range_changed = true;
			if( ! ( ++rt_inner_rng < _inner.count ) ) {
				rt_inner_rng = 0;
				*offo += _outer.range[ rt_outer_rng ].delta;
				if( ! ( *offo < rt_outer_end ) ) {
					if( ++rt_outer_rng < _outer.count )
						outer_range_changed = true;
					else
						return false; // the ONLY case we return false!
				}
			}
		}

		// Above identified the "damage"; now fix it all up.

		if( outer_range_changed ) {
			RANGE_T *pO
				= _outer.range + rt_outer_rng /* the new index */;
			rt_outer_end 
				= is_valid( pO->bounds.stop ) 
				? _minimum( pO->bounds.stop, _N )
				: _N;
			*offo = pO->bounds.start;
			// If outer range changed, inner MUST have too..
			if( ! _begin_inner( *offo, offi ) ) {
#ifdef _UNITTEST_ITER_
				fprintf( stdout, "restart...\n" );
#endif
				goto _entry; // Inner range is empty! Next?...
			}
		} else
		if( inner_range_changed ) {
			if( ! _begin_inner( *offo, offi ) ) {
#ifdef _UNITTEST_ITER_
				fprintf( stdout, "restart...\n" );
#endif
				goto _entry; // Inner range is empty! Next?...
			}
		}

	} while( *offi == *offo /* NEVER allowed */ );

	return true;
}


#ifdef _UNITTEST_ITER_
int main( int argc, char *argv[] ) {

	if( argc > 3 ) {

		const int N
			= atoi( argv[1] );
		const int parse_error
			= iter_parse( argv[2], argv[3] );

		if( ! parse_error ) {
			int o, i;
#ifdef _DEBUG
			iter_dump( "testing: ", stdout );
#endif
			if( iter_begin( N, &o, &i ) ) {
				do {
					printf( "%d\t%d\n", o, i );
					if( o < 0 || i < 0 ) {
						abort();
					}
				} while( iter_next( &o, &i ) );
			} else
				fprintf( stderr, "failed: iter_begin\n" );
		} else
			fprintf( stderr, "failed: iter_parse: %d\n", parse_error );
	} else
		fprintf( stderr, "%s <N> <outer> <inner>\n", argv[0] );

	return EXIT_SUCCESS;
}
#endif

