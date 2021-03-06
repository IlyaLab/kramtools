
/**
  * THE BIG PICTURE
  *
  * This file:
  * 1) (libmtm.a, actually) loads and parses the input matrix
  * 2) iterates through a sequence of feature pairs specified in one of
  *    a variety of ways
  * 3) for each pair it calls a high-level analysis method that selects
  *    and executes an appropriate covariate analysis.
  * 4) routes the analysis results to either
  *    a) immediate output using a specified formatter or
  *    b) a cache to be post-processed (for FDR control) and subsequently
  *       output.
  *
  * Implementation notes:
  *   1) If row labels are present then, by default, they are parsed to
  *      infer the statistical class of the row. Moreover, the parser is
  *      one for "TCGA" format by default. This is so that Sheila et al.'s
  *      scripts need not pass an option for something they consider
  *      "standard." In order to minimize the danger of mis-interpreted
  *      row labels a conservative definition is used (by the
  *      mtm_sclass_by_prefix function). The first two characters of the row
  *      label must match /[BCDFNO][[:punct:]]/.
  *      Otherwise, mtm_sclass_by_prefix returns MTM_STATCLASS_UNKNOWN,
  *      and the class is inferred from the data content of the row.
  * Audit:
  *   1) All error exits return -1. When used as webservice proper status
  *      returns is the least we can do...
  *   2) stdout is ONLY used for data output (at verbosity 0)
  *   3) Analysis failures (because of degeneracy) still return something
  *      the web service can deal with. Striving for data-as-error,
  *      instead of a separate error channel.
  *   4) Insure mtm_resort is called before named rows are used.
  *   5) Every statistical test must fully initialize 1st 4 members of
  *      struct Statistic
  *   6) Casts, ALL casts!
  */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <assert.h>
#include <ctype.h>
#include <err.h>
#include <alloca.h>

#include <gsl/gsl_errno.h>

#include "mtmatrix.h"
#include "mtheader.h"
#include "mtsclass.h"
#include "mterror.h"
#include "featpair.h"
#include "stattest.h"
#include "analysis.h"
#include "varfmt.h"
#include "fixfmt.h"
#include "limits.h"
#include "version.h"

#ifdef HAVE_LUA
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#endif

/***************************************************************************
 * externs
 */

extern int mtm_sclass_by_prefix( const char *token );

extern int get_base10_ints( FILE *fp, int *index, int n );

/***************************************************************************
 * Globals & statics
 */

static const char *MAGIC_SUFFIX         = "-www";
static const char *MAGIC_FORMAT_ID_STD  = "std";
static const char *MAGIC_FORMAT_ID_TCGA = "tcga";

static const char *TYPE_PARSER_INFER    = "auto";

const char *AUTHOR_EMAIL = "rkramer@systemsbiology.org";

/**
  * False Discovery Rate control state.
  */

static double arg_q_value = 0.0;
#define USE_FDR_CONTROL (arg_q_value > 0.0)

/**
  * This is used simply to preclude bloating the FDR cache with results
  * that can't possibly be relevant to the final BH-calculated p-value
  * threshold for FDR.
  * I'm setting this high for safety since it's just an optimization, but
  * it could probably be MUCH smaller.
  */
static double opt_fdr_cache_threshold = 0.5;


// No default on opt_script because looking for a "default.lua" script
// or any *.lua file invites all sorts of confusion with the defaults
// and precedence of other row selection methods.
static const char *opt_script          = NULL;
static bool        opt_header          = true;
static bool        opt_row_labels      = true;
static const char *opt_type_parser     = NULL;
static const char *opt_preproc_matrix  = NULL; // ...or optarg
static       char *opt_single_pair     = NULL; // non-const because it's split

static const char *opt_pairlist_source = NULL;

static const char *NO_ROW_LABELS       = "matrix has no row labels";
static bool        opt_by_name         = false;
#ifdef HAVE_LUA
static const char *DEFAULT_COROUTINE   = "pair_generator";
#endif
static const char *opt_coroutine       = NULL;
static bool        opt_dry_run         = false;
static const char *opt_format          = "tcga";
static bool     opt_warnings_are_fatal = false;

/**
  * Primary output and critical error messages.
  */
#define V_ESSENTIAL (1)

/**
  * Non-critical information that may, nonetheless, be helpful.
  */
#define V_WARNINGS  (2)

/**
  * Purely informational.
  */
#define V_INFO      (3)

static int opt_verbosity = V_ESSENTIAL;

#ifdef _DEBUG
static bool  dbg_silent               = false;
#endif

#define GLOBAL
GLOBAL unsigned  arg_min_cell_count   = 5;
GLOBAL unsigned  arg_min_mixb_count   = 1;
GLOBAL unsigned  arg_min_sample_count = 2; // < 2 NEVER makes sense
#undef  GLOBAL

static const char *opt_na_regex       = NULL; // must be initialized in main

static double      opt_p_value        = 1.0;

/**
  * Eventually it is (was) intended to support delegating row label
  * interpretation--in particular inference of a row's statistical class
  * based on its row label--to a user-provided Lua function.
  * Currently, only the default interpreter built into libmtm is actually
  * supported, which interprets row labels according to "ISB conventions."
  */
static MTM_ROW_LABEL_INTERPRETER _interpret_row_label = mtm_sclass_by_prefix;

static unsigned opt_status_mask       = COVAN_E_MASK;

#ifdef HAVE_LUA
static lua_State *_L = NULL;
/**
  * Delegate inference of a row's statistical class to a user-provided
  * function.
  */
int _lua_statclass_inference( const char *label ) {
	fprintf( stderr, "interpret %s\n", label );
	return MTM_STATCLASS_UNKNOWN;
}

static void _freeLua() {
	lua_close( _L );
}

/**
  * <source> may be literal Lua code or a filename reference.
  */
static int _load_script( const char *source, lua_State *state ) {

	struct stat info;
	return (access( source, R_OK ) == 0) && (stat( source, &info ) == 0)
		? luaL_dofile(   state, source )
		: luaL_dostring( state, source );
}

#endif

static FILE *_fp_output = NULL;
static void (*_emit)( EMITTER_SIG ) = format_tcga;

static bool _sigint_received = false;

/**
  * Records the number of pairwise tests COMPLETED WITHOUT ERROR
  * but not emitted (because they didn't pass the p-value threshold
  * for emission.
  *
  * This is NOT used when FDR control is active in pairwise; rather it is
  * for use in FDR control calculations *outside* the context of pairwise.
  *
  * Note that the total number of tests *attempted* is *not* separately
  * counted because it is assumed that that is known to the caller who
  * set up the run, after all.
  */
static unsigned _insignificant = 0;
static unsigned _untested      = 0;

/**
  * The binary matrix is accessed at runtime through this variable.
  * Notice, in particular, that the analysis code (analysis.cpp)
  * only has access to pairs of rows, NOT the whole matrix.
  * The matrix' structure and implementation is strictly encapsulated
  * within this file, and code in here dolls out just what analysis
  * requires: row pairs.
  */
static struct mtm_matrix _matrix;

static void _freeMatrix( void ) {
	_matrix.destroy( &_matrix );
}

static void _interrupt( int n ) {
	_sigint_received = true;
}

/***************************************************************************
 * Pipeline
 * A) explicit pairs
 *    1. row pair -> analysis -> filtering -> output
 * B) all-pairs or Lua generated pair lists
 *    1. row pair -> analysis -> filtering -> output
 *    2. row pair -> analysis -> fdr aggregation -> re-processing ->output
 */

/**
  * Analysis results on feature pairs are passed EITHER to:
  * 1. an emitter for immediate output (single pair only)
  * 2. a filter (and then possibly to an emitter)
  * 3. an aggregator (to be filtered and stored for later emission)
  *
  * All three of these methods must have the same signature.
  */
#define ANALYSIS_FN_SIG const struct feature_pair *pair

typedef void (*ANALYSIS_FN)( ANALYSIS_FN_SIG );

/***************************************************************************
 * Encapsulates all the decision making regarding actual emission of
 * results.
 */
static void _filter( ANALYSIS_FN_SIG ) {

	struct CovariateAnalysis covan;
	memset( &covan, 0, sizeof(covan) );
	covan_exec( pair, &covan );

	// One last thing to check before filtering to insure corner cases
	// don't fall through the following conditionals....

	if( ! ( isfinite( covan.result.probability ) && fpclassify( covan.result.probability ) != FP_SUBNORMAL ) ) {
		covan.result.probability = 1.0;
		covan.status             = COVAN_E_MATH;
		// The assumption is that if a NaN shows up in the result, it was
		// triggered by an un"pre"detected degeneracy, and so the pair was
		// in fact untestable (how it will be interpreted below).
	}

#ifdef _DEBUG
	if( ! dbg_silent ) {
#endif

		if( ( covan.status & opt_status_mask ) == 0 ) {
			if( covan.result.probability <= opt_p_value )
				_emit( pair, &covan, _fp_output );
			else
				_insignificant += 1;
		} else
			_untested += 1;

#ifdef _DEBUG
	}	
#endif
}

static ANALYSIS_FN _analyze = _filter;

static void _error_handler(const char * reason,
                        const char * file,
                        int line,
                        int gsl_errno) {
	fprintf( stderr,
		"#GSL error(%d): %s\n"
		"#GSL error  at: %s:%d\n",
		gsl_errno, reason, file, line );
}


static FILE *_fdr_cache_fp = NULL;

/***************************************************************************
  * FDR processing
  * This currently involves some redundant computation in the interest of
  * simplicity. Since the whole point of FDR is to limit output, though, the
  * actual runtime cost should be bearabe. Specifically...
  *
  * 1. Two passes are made over the (selected) pairs.
  * 2. The offsets and p-value of the result of the first pass are cached
  *    in a tmp file.
  * 3. After completion, the p-value appropriate for the given q-value is
  *    determined and results from the 1st pass with sufficiently low
  *    p-values are recalculated and, this time, fully emitted.
  *
  * TODO: Optimization: Pairs with very high p-values (i.e. *clearly*
  * uninteresting pairs) can be omitted from the cache even in the first
  * pass and merely counted as tests since their recomputation will *almost*
  * certainly, depending ultimately on the threshold, not be required.
  */

struct FDRCacheRecord {
	double p;
	unsigned a,b;
} __attribute__((packed));
typedef struct FDRCacheRecord FDRCacheRecord_t;
typedef const FDRCacheRecord_t FDRCACHERECORD_T;

static int _cmp_fdr_cache_records( const void *pvl, const void *pvr ) {

	FDRCACHERECORD_T *l = (FDRCACHERECORD_T*)pvl;
	FDRCACHERECORD_T *r = (FDRCACHERECORD_T*)pvr;

	if( l->p == r->p )
		return  0;
	else
		return (l->p < r->p) ? -1 : +1;
}


static int _fdr_uncached_count = 0;


/**
  * Analyze the pair and cache just the offsets and p-value of the result
  * for possible recalculation during post-processing--after FDR control
  * has calculated an appropriate p-value threshold from the q-value.
  */
static void _fdr_cache( ANALYSIS_FN_SIG ) {

	struct CovariateAnalysis covan;
	memset( &covan, 0, sizeof(covan) );

	covan_exec( pair, &covan );

	// Failed tests (for reasons of one kind of degeneracy or another)
	// do not contribute to the calculation of the p-value threshold.

	if( covan.status == 0
		&& isfinite( covan.result.probability ) ) {

		// ...then, whatever the p-value, the test was at least
		// successfully *executed*.

		if( covan.result.probability <= opt_fdr_cache_threshold ) {
			struct FDRCacheRecord rec = {
				.p = covan.result.probability,
				.a = pair->l.offset,
				.b = pair->r.offset
			};
			fwrite( &rec, sizeof(rec), 1, _fdr_cache_fp );
		} else
			_fdr_uncached_count += 1;
	}
}


/**
  * This implements the Benjamini-Hochberg algorithm as described on
  * page 49 of "Large-Scale Inference", Bradley Efron, Cambridge.
  * "...for a fixed value of q in (0,1), let i_max be the largest
  *  index for which
  *                     p_(i) <= (i/N)q
  *
  *  ...and reject H_{0(i)}, the null hypothesis corresponding to
  *  p_(i), if
  *                         i <= i_max,
  *
  *  ...accepting H_{0(i)} otherwise."
  *
  * [...And i is 1-based in this notation!]
  */
static void _fdr_postprocess( FILE *cache, double Q, FILE *final_output, bool minimal_output ) {

	const unsigned CACHED_COUNT
		= ftell( cache ) / sizeof(struct FDRCacheRecord);
	const unsigned TESTED_COUNT
		= CACHED_COUNT
		+ _fdr_uncached_count;

	struct FDRCacheRecord *prec, *sortbuf
		= calloc( CACHED_COUNT, sizeof(struct FDRCacheRecord) );

	const double RATIO
		= Q/TESTED_COUNT;
	int i = 0;

	// Load and sort the cached test records...

	rewind( cache );
	fread( sortbuf, sizeof(struct FDRCacheRecord), CACHED_COUNT, cache );
	qsort( sortbuf, CACHED_COUNT, sizeof(struct FDRCacheRecord), _cmp_fdr_cache_records );

	// ...and recompute the full statistics of all earlier tests that
	// pass the now-established p-value threshold.

	prec = sortbuf;
	if( minimal_output ) {

		while( prec->p <= (i+1)*RATIO && i < CACHED_COUNT) 
			fprintf( final_output, "%d\t%d\t%.3e\n", prec->a, prec->b, prec->p );

	} else {

		while( prec->p <= (i+1)*RATIO && i < CACHED_COUNT) {
		
			struct feature_pair fpair;
			struct CovariateAnalysis covan;
			memset( &covan, 0, sizeof(covan) );

			fpair.l.offset = prec->a;
			fpair.r.offset = prec->b;
			fetch_by_offset( &_matrix, &fpair );

			covan_exec( &fpair, &covan );

			// At this point emission is unconditional; FDR control has
			// already filtered all that will be filtered...

			_emit( &fpair, &covan, final_output );

			if( _sigint_received ) {
				time_t now = time(NULL);
				fprintf( stderr, "# FDR postprocess interrupted @ %s", ctime(&now) );
				break;
			}

			prec += 1;
			i    += 1;
		}
	}

	// However, the preceding loop exited

	if( opt_verbosity >= V_WARNINGS ) {
		if( i > 0 )
			fprintf( final_output, "# max p-value %.3f\n", sortbuf[i-1].p );
		else
			fprintf( final_output, "# no values passed FDR control\n" );
	}
	if( sortbuf )
		free( sortbuf );
}


/***************************************************************************
  * Row selection iterators
  * Each of these 5 methods takes arguments specific to the
  * source of input pairs and uses file static state for output.
  * All but the first delegate to the _analyze function which may be a
  * filter on immediate-mode output or a caching function for FDR.
  */

// BEGIN:RSI

static int /*ANCP*/ _analyze_cross_product(
		const struct mtm_matrix_header *hdr, FILE *fp[] ) {

	bool completed = true;
	struct feature_pair fpair;
	struct mtm_row *rrid;

	/**
	  * Location of left data won't change: same buffer, same offset for
	  * duration of iteration.
	  */
	fpair.l.data
		= calloc( _matrix.columns, sizeof(mtm_int_t) );

	/**
	  * TODO: I actually could enumerate the disk-resident matrix'
	  * row names just as I do the descriptors. In general, the
	  * pervasive assumption throughout this source that the input
	  * is exactly one matrix needs to be revisited.
	  */
	fpair.l.name = NULL;

	if( fpair.l.data == NULL )
		return -1;

	assert( ! _matrix.lexigraphic_order /* should be row order */ );

	for(fpair.l.offset = 0;
		fpair.l.offset < hdr->rows;
		fpair.l.offset++ ) {

		/**
		  * Read the "left" feature's data and descriptor
		  */

		if( fread( (void*)fpair.l.data, sizeof(mtm_int_t), hdr->columns,  fp[0] ) != hdr->columns )
			break;
		if( fread( &fpair.l.desc, sizeof(struct mtm_descriptor), 1, fp[1] ) != 1 )
			break;

		/**
		  * RAM-resident matrix is *fully* reset for each row of disk-
		  * resident matrix...
		  */

		fpair.r.data = _matrix.data;
		rrid         = _matrix.row_map; // may be NULL

		for(fpair.r.offset = 0;
			fpair.r.offset < _matrix.rows;
			fpair.r.offset++ ) {

			fpair.r.name = rrid ? rrid->string : "";
			fpair.r.desc = _matrix.desc[ fpair.r.offset ];

			_analyze( &fpair );

			if( _sigint_received ) {
				time_t now = time(NULL);
				fprintf( stderr, "# main analysis loop interrupted @ %s", ctime(&now) );
				completed = false;
				break;
			}

			fpair.r.data += _matrix.columns;
			if( rrid ) rrid += 1;

		} // inner for
	}

	if( ferror( fp[0] ) || ferror( fp[1] ) ) {
		warn( "reading row %d of preprocessed matrix", fpair.l.offset );
		completed = false;
	}

	if( fpair.l.data )
		free( (void*)fpair.l.data );

	return completed ? 0 : -1;
}


static bool _is_integer( const char *pc ) {
	while( *pc ) if( ! isdigit(*pc++) ) return false;
	return true;
}


static int _analyze_single_pair( const char *csv, const bool HAVE_ROW_LABELS ) {

	static const char *MISSING_MSG
		= "you specified row %s for a matrix without row names.\n";

	int econd;
	char *right, *left = alloca( strlen(csv)+1 );
	struct CovariateAnalysis covan;
	struct feature_pair pair;

	memset( &pair,  0, sizeof(pair) );
	memset( &covan, 0, sizeof(covan) );

	// Split the string in two...

	strcpy( left, csv );
	right = strchr( left, ',' );
	if( NULL == right )
		errx( -1, "missing comma separator in feature pair \"%s\"", left );
	*right++ = '\0';

	// Lookup each part.(They need not be the same format.)

	if( _is_integer( left ) ) {
		pair.l.offset = atoi( left );
		mtm_resort_rowmap( &_matrix, MTM_RESORT_BYROWOFFSET );
		econd = mtm_fetch_by_offset( &_matrix, &(pair.l) );
	} else
	if( HAVE_ROW_LABELS ) {
		pair.l.name   = left;
		if( mtm_resort_rowmap( &_matrix, MTM_RESORT_LEXIGRAPHIC ) )
			errx( -1, NO_ROW_LABELS );
		econd = mtm_fetch_by_name( &_matrix, &(pair.l) );
	} else
		errx( -1, MISSING_MSG, left );

	if( _is_integer( right ) ) {
		pair.r.offset = atoi( right );
		mtm_resort_rowmap( &_matrix, MTM_RESORT_BYROWOFFSET );
		econd = mtm_fetch_by_offset( &_matrix, &(pair.r) );
	} else
	if( HAVE_ROW_LABELS ) {
		pair.r.name = right;
		if( mtm_resort_rowmap( &_matrix, MTM_RESORT_LEXIGRAPHIC ) )
			errx( -1, NO_ROW_LABELS );
		econd = mtm_fetch_by_name( &_matrix, &(pair.r) );
	} else
		errx( -1, MISSING_MSG, right );

	covan_exec( &pair, &covan );

	_emit( &pair, &covan, _fp_output );

	return 0;
}


static int /*ANAM*/ _analyze_named_pair_list( FILE *fp ) {

	struct feature_pair fpair;

	size_t blen = 0;
	ssize_t llen;
	char *right, *left = NULL;

	while( ( llen = getline( &left, &blen, fp ) ) > 0 ) {

		// Trim off the newline char.

		if( left[llen-1] == '\n' )
			left[--llen] = '\0'; // Strip the NL

		// Parse the two names out of the line.

		right = strchr( left, '\t' );
		if( right )
			*right++ = '\0';
		else {
			fprintf( stderr, "error: no tab found in '%s'.\n", left );
			if( opt_warnings_are_fatal )
				break;
			else
				continue; // no reason we -can't- continue
		}

		fpair.l.name  = left;
		fpair.r.name = right;

		if( fetch_by_name( &_matrix, &fpair ) ) {

			warnx( "error: one or both of...\n"
				"\t1) %s\n"
				"\t2) %s\n"
				"\t...not found.\n",
				fpair.l.name,
				fpair.r.name );
			if( opt_warnings_are_fatal )
				break;
			else
				continue; // no reason we -can't- continue
		} else
			_analyze( &fpair );
	}

	if( left )
		free( left );

	return 0;
}


static int /*ANUM*/ _analyze_pair_list( FILE *fp ) {

	struct feature_pair fpair;

	int arr[2];

	while( get_base10_ints( fp, arr, 2 ) == 2 ) {

		fpair.l.offset = arr[0];
		fpair.r.offset = arr[1];

		if( fetch_by_offset( &_matrix, &fpair ) ) {
			warnx( "error: one of row indices (%d,%d) not in [0,%d)\n"
				"\tjust before byte offset %ld in the stream.\n"
				"\tAborting...\n",
				fpair.l.offset,
				fpair.r.offset,
				_matrix.rows,
				ftell( fp ) );
			if( opt_warnings_are_fatal )
				break;
			else	
				continue; // no reason we -can't- continue
		} else
			_analyze( &fpair );
	}
	return 0;
}


#ifdef HAVE_LUA
static int /*ALUA*/ _analyze_generated_pair_list( lua_State *state ) {

	struct feature_pair fpair;
	int isnum, lua_status;

	lua_getglobal( _L, opt_coroutine );

	assert( ! lua_isnil( _L, -1 ) /* because it was checked early */ );

	do {

		lua_pushnumber( _L, _matrix.rows );
		lua_status = lua_resume( _L, NULL, 1 );

		if( lua_status == LUA_YIELD ) {

			fpair.r.offset = lua_tonumberx( _L, -1, &isnum );
			fpair.l.offset = lua_tonumberx( _L, -2, &isnum );
			lua_pop( _L, 2 );

			if( fetch_by_offset( &_matrix, &fpair ) ) {
				warnx( "one or both of %s-generated row indices (%d,%d) not in [0,%d)\n",
					opt_coroutine,
					fpair.l.offset,
					fpair.r.offset,
					_matrix.rows );
				if( opt_warnings_are_fatal )
					break;
				else	
					continue; // no reason we -can't- continue
			} else
				_analyze( &fpair );

		} else
		if( lua_status == LUA_OK )
			break;
		else { // some sort of error occurred.
			fputs( lua_tostring( _L, -1 ), stderr );
		}

		if( _sigint_received ) {
			time_t now = time(NULL);
			fprintf( stderr, "analysis loop interrupted @ %s", ctime(&now) );
			break;
		}

	} while( lua_status == LUA_YIELD );


	return 0;
}
#endif

/**
  * This clause serves the primary use case motivating this
  * application: FAST, EXHAUSTIVE (n-choose-2) pairwise analysis.
  * As a result, the coding of this iteration schema is quite
  * different from the others. In particular:
  * 1. Since I -know- the order of feature pair evaluation I can
  *    preclude lots of useless work by entirely skipping outer loop
  * features with univariate degeneracy. This isn't feasible in
  * the other iteration sc
  * I'm dispensing with "pretty" and doing the
  * pointer arithmetic locally to eek out maximum possible speed.
  * In particular, I set up lots of ptrs below to obviate
  * redundant pointer arithmetic--only do fixed additions.
  */
static int /*AALL*/ _analyze_all_pairs( void ) {

	bool completed = true;
	struct feature_pair fpair;

	struct mtm_row *lrid, *rrid;

	assert( ! _matrix.lexigraphic_order /* should be row order */ );

	fpair.l.data = _matrix.data;
	lrid         = _matrix.row_map; // may be NULL

	for(fpair.l.offset = 0;
		fpair.l.offset < _matrix.rows;
		fpair.l.offset++ ) {

		fpair.l.name = lrid ? lrid->string : "";
		fpair.l.desc = _matrix.desc[ fpair.l.offset ];

		fpair.r.data = fpair.l.data + _matrix.columns;
		rrid = lrid ? lrid + 1 : NULL;

		for(fpair.r.offset = fpair.l.offset+1;
			fpair.r.offset < _matrix.rows;
			fpair.r.offset++ ) {

			fpair.r.name = rrid ? rrid->string : "";
			fpair.r.desc = _matrix.desc[ fpair.r.offset ];

			_analyze( &fpair );

			if( _sigint_received ) {
				time_t now = time(NULL);
				fprintf( stderr, "# main analysis loop interrupted @ %s", ctime(&now) );
				completed = false;
				break;
			}

			fpair.r.data += _matrix.columns;
			if( rrid ) rrid += 1;

		} // inner for

		fpair.l.data += _matrix.columns;
		if( lrid ) lrid += 1;
	}
	return completed ? 0 : -1;
}

// END:RSI

/**
  * Initializations for which static initialization can't/shouldn't 
  * be relied upon. Common to executable and Python extension.
  */
static void _jit_initialization( void ) {
	opt_na_regex = mtm_default_NA_regex;
	memset( &_matrix,  0, sizeof(struct mtm_matrix) );
	gsl_set_error_handler( _error_handler );
}


/***************************************************************************
  * Online help
  */

static const char *_YN( bool y ) {
	return y ? "yes" : "no";
}

#define USAGE_SHORT false
#define USAGE_LONG  true

static void _print_usage( const char *exename, FILE *fp, bool exhaustive ) {

	extern const char *USAGE_UNABRIDGED;
	extern const char *USAGE_ABRIDGED;

	char debug_state[64];
	debug_state[0] = 0;
#ifdef _DEBUG
	sprintf( debug_state, "(DEBUG) Compiled %s %s", __DATE__,__TIME__ );
#endif

	if( exhaustive )
		fprintf( fp, USAGE_UNABRIDGED,
			exename, VER_MAJOR, VER_MINOR, VER_PATCH, VER_TAG, debug_state,
			exename,
			exename,
			TYPE_PARSER_INFER,
			opt_na_regex,
#ifdef HAVE_LUA
			DEFAULT_COROUTINE,
#endif
			arg_min_cell_count,
			arg_min_mixb_count,
			arg_min_sample_count,

			opt_p_value,
			opt_status_mask,
			opt_format,
			MAGIC_FORMAT_ID_STD, MAGIC_FORMAT_ID_TCGA,
			_YN(opt_warnings_are_fatal),
			opt_verbosity,
			MAGIC_SUFFIX,
			MAX_CATEGORY_COUNT,
			MTM_MAX_MISSING_VALUES,
			AUTHOR_EMAIL );
	else
		fprintf( fp, USAGE_ABRIDGED,
			exename, VER_MAJOR, VER_MINOR, VER_PATCH, VER_TAG, debug_state,
			exename,
			opt_p_value,
		  	AUTHOR_EMAIL );
}


int main( int argc, char *argv[] ) {

	int exit_status = EXIT_SUCCESS;

	const char *i_file = NULL;
	const char *o_file = NULL;
	FILE *fp           = NULL;

	if( argc < 2 ) { // absolute minimum args: <executable name> <input matrix>
		_print_usage( argv[0], stdout, USAGE_SHORT );
		exit( EXIT_SUCCESS );
	}

	_jit_initialization();

	do {

		static const char *CHAR_OPTIONS
#ifdef HAVE_LUA
			= "s:hrt:N:C:P:n:x:c:DM:p:f:q:v:?X";
#else
			= "hrt:N:C:P:n:x:DM:p:f:q:v:?X";
#endif

		static struct option LONG_OPTIONS[] = {
#ifdef HAVE_LUA
			{"script",        required_argument,  0,'s'},
#endif
			{"no-header",     no_argument,        0,'h'},
			{"no-row-labels", no_argument,        0,'r'},
			{"type-parser",   required_argument,  0,'t'},
			{"na-regex",      required_argument,  0,'N'},

			{"crossprod",     required_argument,  0,'C'},
			{"pair",          required_argument,  0,'P'},
			{"by-name",       required_argument,  0,'n'},
			{"by-index",      required_argument,  0,'x'},
#ifdef HAVE_LUA
			{"coroutine",     required_argument,  0,'c'},
#endif
			{"dry-run",       no_argument,        0,'D'},

			{"min-ct-cell",   required_argument,  0, 256 }, // no short equivalents
			{"min-mx-cell",   required_argument,  0, 257 }, // no short equivalents
			{"min-samples",   required_argument,  0,'M'},

			{"p-value",       required_argument,  0,'p'},
			{"format",        required_argument,  0,'f'},
			{"fdr",           required_argument,  0,'q'},
			{"verbosity",     required_argument,  0,'v'},
#ifdef _DEBUG
			{"debug",         required_argument,  0, 258 }, // no short equivalents
#endif
			{"help",          no_argument,        0,'?'},
			{ NULL,           0,                  0, 0 }
		};

		int arg_index = 0;
		const int c
			= getopt_long( argc, argv, CHAR_OPTIONS, LONG_OPTIONS, &arg_index );

		switch (c) {

		case 's': // script
			opt_script      = optarg;
			break;
		case 'h': // no-header
			opt_header      = false;
			break;
		case 'r': // no-row-labels
			opt_row_labels  = false;
			break;
		case 't': // type-parser
			opt_type_parser = optarg;
			break;
		case 'N': // na-regex
			opt_na_regex    = optarg;
			break;
		case 'C': // cross-product of matrices
			opt_preproc_matrix = optarg;
			break;
		case 'P': // pair
			opt_single_pair = optarg;
			break;
		case 'n': // by-name
			opt_pairlist_source = optarg;
			opt_by_name         = true;
			break;
		case 'x': // by-index
			opt_pairlist_source = optarg;
			opt_by_name         = false;
			break;
		case 'c': // coroutine
			opt_coroutine       = optarg;
			break;
		case 'D': // dry-run
			opt_dry_run         = true;
			break;

		////////////////////////////////////////////////////////////////////
		case 256: // ...because I haven't defined a short form for this
			arg_min_cell_count = atoi( optarg );
			break;

		case 257: // ...because I haven't defined a short form for this
			arg_min_mixb_count = atoi( optarg );
			break;

		case 'M':
			arg_min_sample_count = atoi( optarg );
			if( arg_min_sample_count < 2 ) {
				warnx( "Seriously...%d samples is acceptable?\n"
					"I don't think so... ;)\n",
					arg_min_sample_count );
				exit( EXIT_FAILURE );
			}
			break;
		////////////////////////////////////////////////////////////////////

		case 'p':
			opt_p_value = atof( optarg );
			if( ! ( 0 < opt_p_value ) ) {
				warnx( "specified p-value %.3f will preclude all output.\n",
					opt_p_value );
				abort();
			} else
			if( ! ( opt_p_value < 1.0 ) ) {
				warnx( "p-value %.3f will filter nothing.\n"
					"\tIs this really what you want?\n",
					opt_p_value );
			}
			break;

		case 'J': // JSON format
		case 'f': // tabular format
			// Check for magic-value strings first
			if( strcmp( MAGIC_FORMAT_ID_STD, optarg ) == 0 )
				_emit = format_standard;
			else
			if( strcmp( MAGIC_FORMAT_ID_TCGA, optarg ) == 0 )
				_emit = format_tcga;
			else {
				const char *specifier
					= emit_config( optarg, c=='J' ? FORMAT_JSON : FORMAT_TABULAR );
				if( specifier ) {
					errx( -1, "invalid specifier \"%s\"", specifier );
				}
				_emit = emit_exec;
			}
			break;

		case 'q':
			arg_q_value = atof( optarg );
			_analyze = _fdr_cache;
			break;

		case 'v': // verbosity
			opt_verbosity = atoi( optarg );
			break;

		case '?': // help
			_print_usage( argv[0], stdout, USAGE_SHORT );
			exit( EXIT_SUCCESS );
			break;

		case 'X': // help
			_print_usage( argv[0], stdout, USAGE_LONG );
			exit( EXIT_SUCCESS );
			break;

#ifdef _DEBUG
		case 258:
			if( strchr( optarg, 'X' ) != NULL ) {
				// Turn OFF all filtering to turn on "exhaustive"
				// output mode.
				opt_status_mask = 0;
				opt_p_value     = 1.0;
			}
			dbg_silent     = strchr( optarg, 'S' ) != NULL;
			break;
#endif
		case -1: // ...signals no more options.
			break;
		default:
			printf ("error: unknown option: %c\n", c );
			exit( EXIT_FAILURE );
		}
		if( -1 == c ) break;

	} while( true );

#ifdef HAVE_LUA

	if( opt_script ) {

		/**
		  * Create a Lua state machine and load/run the Lua script from file
		  * or command line before anything else.
		  */

		_L = luaL_newstate();
		if( _L ) {
			atexit( _freeLua ); // ...so lua_close needed NOWHERE else.
			luaL_openlibs( _L );
			if( _load_script( opt_script, _L ) != LUA_OK ) {
				errx( -1, "in \"%s\": %s",
					opt_script, lua_tostring( _L,-1) );
			}
		} else
			errx( -1, "failed creating Lua statespace" );

		/**
		  * Fail early!
		  * If the user provided a coroutine name, check for its presence.
		  * Otherwise, see if the default coroutine is present.
		  * Otherwise, Lua is not being used to generate pairs.
		  */

		if( opt_coroutine ) {

			lua_getglobal( _L, opt_coroutine );
			if( lua_isnil( _L, -1 ) ) {
				errx( -1, "pair generator coroutine \"%s\" not defined in \"%s\"",
					opt_coroutine, opt_script );
			}

		} else {

			lua_getglobal( _L, DEFAULT_COROUTINE );

			if( ! lua_isnil( _L, -1 ) )
				opt_coroutine = DEFAULT_COROUTINE;

			// Non-existence of DEFAULT_COROUTINE is not an error;
			// we assume user does NOT intend Lua to generate the pairs.
		}

		lua_pop( _L, 1 ); // clean the stack

		/**
		  * Similarly verify NOW that the a specified type parser is present.
		  */
		if( opt_type_parser ) {
			lua_getglobal( _L, opt_type_parser );
			if( lua_isnil( _L, -1 ) ) {
				errx( -1, "pair generator coroutine \"%s\" not defined in \"%s\"",
					opt_type_parser, opt_script );
			}
			lua_pop( _L, 1 ); // clean the stack
		}

		/**
		  * Load Lua script and execute a dry-run if applicable.
		  * This is in advance of other argument validation since other
		  * arguments will be ignored...
		  */

		if( opt_dry_run && opt_coroutine != NULL ) {

			const int COUNT
				= getenv("DRY_RUN_COUNT")
				? atoi( getenv("DRY_RUN_COUNT") )
				: 8;

			int isnum, lua_status = LUA_YIELD;

			lua_getglobal( _L, opt_coroutine );
			if( lua_isnil( _L, -1 ) )
				errx( -1, "%s not defined (in Lua's global namespace)",
					opt_coroutine );

			while( lua_status == LUA_YIELD ) {
				lua_pushnumber( _L, COUNT );
				lua_status = lua_resume( _L, NULL, 1 );
				if( lua_status <= LUA_YIELD /* OK == 0, YIELD == 1*/ ) {
					// Only output coroutine.yield'ed values...
					if( lua_status == LUA_YIELD ) {
						const int b = lua_tonumberx( _L, -1, &isnum );
						const int a = lua_tonumberx( _L, -2, &isnum );
						lua_pop( _L, 2 );
						fprintf( stdout, "%d %d\n", a, b );
					}
				} else
					fputs( lua_tostring( _L, -1 ), stderr );
			}

			exit( EXIT_SUCCESS );
		}
	} else
	if( opt_coroutine ) {
		errx( -1, "coroutine specified (\"%s\") but no script", opt_coroutine );
	}

#endif

	/**
	  * Catch simple argument inconsistencies and fail early!
	  */

	if( opt_single_pair ) {
		if( USE_FDR_CONTROL ) {
			warnx( "FDR is senseless on a single pair.\n" );
			if( opt_warnings_are_fatal )
				abort();
			else
				arg_q_value = 0.0;
		}
	}

	/**
	  * The last two positional arguments are expected to be filenames
	  * ... <filename1>
	  */

	static const char *NAME_STDIN  = "stdin";
	static const char *NAME_STDOUT = "stdout";

	switch( argc - optind ) {

	case 0: // input MUST be stdin, output stdout
		i_file = NAME_STDIN;
		o_file = NAME_STDOUT;
		break;

	case 1:
		if( access( argv[ optind ], R_OK ) == 0 ) {
			i_file = argv[ optind++ ];
			o_file = NAME_STDOUT;
		} else {
			i_file = NAME_STDIN;
			o_file = argv[ optind++ ];
		}
		break;

	case 2:
	default:
		i_file = argv[ optind++ ];
		o_file = argv[ optind++ ];
		if( (argc-optind) > 0 && opt_verbosity >= V_ESSENTIAL ) {
			// This behavior of ignoring "extra" arguments is implemented
			// ONLY so that this executable plays nicely with a job control
			// system which uses extra command line args NOT intended for
			// the executable. TODO: Revisit this!
			fprintf( stderr,
				"warning: ignoring %d trailing positional arguments.\n", argc-optind );
			while( optind < argc )
				fprintf( stderr, "\t\"%s\"\n", argv[optind++] );
			if( opt_warnings_are_fatal )
				exit( EXIT_FAILURE );
		}
	}

	if( opt_pairlist_source != NULL
			&& strcmp( opt_pairlist_source, i_file ) == 0 ) {
		errx( -1, "error: stdin specified (or implied) for both pair list and the input matrix\n" );
	}

	/**
	  * Emit intentions...
	  */

	if( opt_verbosity >= V_INFO ) {

#define MAXLEN_FS 75

		char feature_selection[MAXLEN_FS+1];
		feature_selection[MAXLEN_FS] = 0; // insure NUL termination

		/**
		  * Following conditional cascade must exactly match the larger
		  * one below in order to accurately report feature selection method.
		  */
		if( opt_single_pair ) {
			strncpy( feature_selection, opt_single_pair, MAXLEN_FS );
		} else
		if( opt_pairlist_source ) {
			snprintf( feature_selection,
				MAXLEN_FS,
				"by %s in %s",
				opt_by_name ? "name" : "offset", opt_pairlist_source );
		} else {
#ifdef HAVE_LUA
			if( opt_coroutine )
				snprintf( feature_selection,
					MAXLEN_FS,
					"%s in %s",
					opt_coroutine, opt_script /* may be literal source */ );
			else
#endif
				strncpy( feature_selection, "all-pairs", MAXLEN_FS );
		}

		fprintf( stderr,
			"       input: %s\n"
			"      select: %s\n"
			"      output: %s\n"
#ifdef _DEBUG
			"       debug: silent: %s\n"
#endif
			,i_file,
			feature_selection,
			o_file
#ifdef _DEBUG
			, _YN( dbg_silent )
#endif
			);
	}

	/**
	  * Load the input matrix.
	  */

	fp = strcmp( i_file, NAME_STDIN )
		? fopen( i_file, "r" )
		: stdin;
	if( fp ) {

		const unsigned int FLAGS
			= ( opt_header ? MTM_MATRIX_HAS_HEADER : 0 )
			| ( opt_row_labels ? MTM_MATRIX_HAS_ROW_NAMES : 0 )
			| ( opt_verbosity & MTM_VERBOSITY_MASK);

		const int econd
			= mtm_parse( fp,
				FLAGS,
				opt_na_regex,
				MAX_CATEGORY_COUNT,
				opt_row_labels ? _interpret_row_label : NULL,
				NULL, // ...since no persistent binary matrix is needed.
				&_matrix );
		fclose( fp );

		if( econd )
			errx( -1, "mtm_parse returned (%d)", econd );
		else
			atexit( _freeMatrix );
	}

	if( opt_dry_run ) { // a second possible
		exit( EXIT_SUCCESS );
	}

	if( SIG_ERR == signal( SIGINT, _interrupt ) ) {
		warn( "failed installing interrupt handler\n"
			"\tCtrl-C will terminated gracelessly\n" );
	}

	/**
	  * Choose and open, if necessary, an output stream.
	  */

	_fp_output
		= (strcmp( o_file, NAME_STDOUT ) == 0 )
		? stdout
		: fopen( o_file, "w" );

	if( NULL == _fp_output ) {
		err( -1, "opening output file \"%s\"", o_file );
	}

	if( covan_init( _matrix.columns ) ) {
		err( -1, "error: covan_init(%d)\n", _matrix.columns );
	} else
		atexit( covan_fini );

	if( opt_verbosity >= V_INFO )
		fprintf( _fp_output, "# %d rows/features X %d columns/samples\n", _matrix.rows, _matrix.columns );

	if( USE_FDR_CONTROL ) {
		_fdr_cache_fp = tmpfile();
		_fdr_uncached_count = 0;
		if( NULL == _fdr_cache_fp )
			err( -1, "creating a temporary file" );
	}

	/**
	  * Here the main decision is made regarding feature selection.
	  * In order of precedence:
	  * 0. cross-product of matrices
	  * 1. single pair
	  * 2. explicit pairs (by name or by offset)
	  * 3. Lua-generated offsets
	  * 4. all-pairs
	  */

	if( opt_preproc_matrix ) {
		FILE *ppm[2];
		ppm[0] = fopen( opt_preproc_matrix, "r" );
		if( ppm[0] ) {

			struct mtm_matrix_header hdr;

			/**
			  * Open a second stream on the file to access the
			  * descriptor table.
			  */
			ppm[1] = fopen( opt_preproc_matrix, "r" );
			
			/**
			  * Read the header and verify compatibility
			  */

			if( mtm_load_header( ppm[0], &hdr ) != MTM_OK )
				err( -1, "failed reading preprocessed matrix' (%s) header",
					opt_preproc_matrix );

			// Verify file is the preprocessed matrix the user thinks it is.

			if( strcmp( hdr.sig, MTM_SIGNATURE ) )
				errx( -1, "%s has wrong signature."
					"Are you sure this is a preprocessed matrix",
					opt_preproc_matrix );

			// Verify equality of columns

			if( _matrix.columns != hdr.columns )
				errx( -1, "processed matrix has %d columns; other has %d",
					hdr.columns, _matrix.columns );

			// Skip past the header to the row data.

			if( fseek( ppm[0], hdr.section[ S_DATA ].offset, SEEK_SET ) )
				err( -1, "failed seeking to start of data" );

			// Skip to the descriptor table in the 2nd stream.

			if( fseek( ppm[1], hdr.section[ S_DESC ].offset, SEEK_SET ) )
				err( -1, "failed seeking to start of data" );

			_analyze_cross_product( &hdr, ppm );

			fclose( ppm[1] );
			fclose( ppm[0] );
		}
	} else
	if( opt_single_pair ) {

		_analyze_single_pair( opt_single_pair, opt_row_labels );

	} else
	if( opt_pairlist_source ) {

		FILE *fp
			= strcmp(opt_pairlist_source,NAME_STDIN)
			? fopen( opt_pairlist_source, "r" )
			: stdin;
	
		if( opt_by_name ) {	
			if( mtm_resort_rowmap( &_matrix, MTM_RESORT_LEXIGRAPHIC ) )
				errx( -1, NO_ROW_LABELS );
		}

		if( fp ) {
			const int econd
				= opt_by_name
				? _analyze_named_pair_list( fp )
				: _analyze_pair_list( fp );
			if( econd )
				warn( "error (%d) analyzing%s pair list",
					econd, opt_by_name ? " named" : "" );
			fclose( fp );
		} else
			warn( "opening \"%s\"", opt_pairlist_source );

	} else {

#ifdef HAVE_LUA
		if( opt_coroutine )
			_analyze_generated_pair_list( _L );
		else
#endif
			_analyze_all_pairs();
	}

	// Post process results if FDR is in effect and the 1st pass was
	// allowed to complete. (Post-processing involves a repetition of
	// analysis FOR A SUBSET of the original input.)

	if( USE_FDR_CONTROL && (! _sigint_received) ) {
		_fdr_postprocess( _fdr_cache_fp, arg_q_value, _fp_output, opt_preproc_matrix != NULL );
	} else
	if( opt_verbosity >= V_ESSENTIAL ) {
		fprintf( _fp_output, 
				"# %d filtered for insignificance\n"
				"# %d filtered for some sort of degeneracy\n", 
				_insignificant,
				_untested );
		// ...which does not apply in FDR control context.
	}

	if( _fdr_cache_fp )
		fclose( _fdr_cache_fp );

	if( _fp_output )
		fclose( _fp_output );

	return exit_status;
}

