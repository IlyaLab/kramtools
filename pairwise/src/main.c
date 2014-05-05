
/**
  * This file serves 3 primary functions:
  * 1) loading and parsing the input matrix (done by libmtm.a) 
  * 2) iterating through a sequence of feature pairs (specified
  *    in a variety of ways) and for each pair calling a high-level 
  *    analysis method that selects an appropriate covariate analysis.
  * 3) routing the results to either
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

#include <gsl/gsl_errno.h>

#include "mtmatrix.h"
#include "mtsclass.h"
#include "featpair.h"
#include "stattest.h"
#include "analysis.h"
#include "varfmt.h"
#include "fixfmt.h"
#include "limits.h"

#ifdef HAVE_LUA
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#endif

/***************************************************************************
 * externs
 */

extern int get_base10_ints( FILE *fp, int *index, int n );
extern int mtm_sclass_by_prefix( const char *token );

/***************************************************************************
 * Globals & statics
 */

static const char *MAGIC_SUFFIX = "-www";
static const char *MAGIC_STD  = "std";
static const char *MAGIC_TCGA = "tcga";

const char *AUTHOR_EMAIL = "rkramer@systemsbiology.org";

FILE *g_fp_output = NULL;
FILE *g_fp_cache = NULL;

size_t g_COLUMNS = 0;

static void (*_emit)( EMITTER_SIG ) = format_tcga;

#define GLOBAL
GLOBAL unsigned  arg_min_cell_count   = 5;
GLOBAL unsigned  arg_min_mixb_count   = 1;
GLOBAL unsigned  arg_min_sample_count = 2; // < 2 NEVER makes sense
GLOBAL bool g_sigint_received          = false; // fdr.d needs access
#undef  GLOBAL

static const char *TYPE_PARSER_INFER   = "auto";
static double arg_q_value              = 0.0;
#define USE_FDR_CONTROL (arg_q_value > 0.0)

// No default on opt_script because looking for a "default.lua" script
// or any *.lua file invites all sorts of confusion with the defaults
// and precedence of other row selection methods.
static const char *opt_script          = NULL;
static bool        opt_header          = true;
static bool        opt_row_labels      = true;
static const char *opt_type_parser     = NULL;
static const char *opt_na_regex        = NULL; // initialized in main

static       char *opt_single_pair     = NULL; // non-const because it's split

static const char *opt_pairlist_source = NULL;
static bool        opt_by_name         = false;
static const char *DEFAULT_COROUTINE   = "pair_generator";
static const char *opt_coroutine       = NULL;
static bool        opt_dry_run         = false;

static double      opt_p_value         = 1.0;
static const char *opt_format          = "tcga";

static bool     opt_warnings_are_fatal = false;

static int         opt_verbosity       = 0;

static bool  opt_running_as_webservice = false;

static MTM_ROW_LABEL_INTERPRETER _interpret_row_label = mtm_sclass_by_prefix;

static unsigned opt_status_mask        = COVAN_E_MASK;

#ifdef _DEBUG
/**
  * Emit all results INCLUDING those with degeneracy.
  */

static bool  dbg_silent     = false;
#endif

#ifdef HAVE_LUA
static lua_State *_L = NULL;
int _delegate_row_label_to_lua( const char *label ) {
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


/**
  * Records the number of each test type that are filtered.
  */
// TODO: static unsigned _filtered[ CovarTypeCount ];

/**
  * The binary matrix is accessed at runtime through this variable.
  * Notice, in particular, that the analysis code (analysis.cpp)
  * only has access to pairs of rows, NOT the whole matrix.
  * The matrix' structure and implementation is strictly encapsulated
  * within this file, and code in here dolls out just what analysis
  * requires: row pairs.
  */
static struct mtm_matrix _matrix;

static void _interrupt( int n ) {
	g_sigint_received = true;
}


static void _freeMatrix() {
	_matrix.destroy( &_matrix );
}


void panic( const char *src, int line ) {
	fprintf( stderr, "panic on %s:%d", src, line );
	abort();
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

	if( ! isfinite( covan.result.probability ) ) { 
		covan.result.probability = 1.0;
		covan.status             = COVAN_E_MATH;
	}

#ifdef _DEBUG
	if( ! dbg_silent ) {
#endif
		if( ( covan.result.probability <= opt_p_value ) 
			&& 
			( ( covan.status & opt_status_mask ) == 0 ) ) {

			_emit( pair, &covan, g_fp_output );

		} else {
			//_filtered[ g_summary.kind ] += 1;
		}

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


/***************************************************************************
  * FDR processing
  * This 
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


/**
  */
static void _fdr_cache( ANALYSIS_FN_SIG ) {

	struct CovariateAnalysis covan;
	memset( &covan, 0, sizeof(covan) );
	covan_exec( pair, &covan );

	if( covan.status == 0 && covan.result.probability <= opt_p_value ) {
		struct FDRCacheRecord rec = {
			.p = covan.result.probability,
			.a = pair->l.offset,
			.b = pair->r.offset
		};
		fwrite( &rec, sizeof(rec), 1, g_fp_cache );
	}
}


static void fdr_postprocess( FILE *cache, double Q ) {

	const long size
		= ftell( cache );
	const unsigned COUNT
		= size / sizeof(struct FDRCacheRecord);
	struct FDRCacheRecord *prec, *sortbuf
		= calloc( size, sizeof(struct FDRCacheRecord) );
	const double RATIO
		= Q/COUNT;
	int n = 0;

	rewind( cache );
		fread( sortbuf, sizeof(struct FDRCacheRecord), COUNT, cache );
	rewind( cache );

	qsort( sortbuf, COUNT, sizeof(struct FDRCacheRecord), _cmp_fdr_cache_records );

	prec = sortbuf;
	while( prec->p <= (n+1)*RATIO ) {

		// TODO:const unsigned *pa = g_matrix_body + g_COLUMNS*prec->a;
		// TODO:const unsigned *pb = g_matrix_body + g_COLUMNS*prec->b;

/*
		struct CovariateAnalysis covan;
		memset( &covan, 0, sizeof(covan) );
		covan_exec( pair, &covan );
*/
		//_emit( pair, &covan, g_fp_output );

		if( g_sigint_received ) {
			time_t now = time(NULL);
			fprintf( stderr, "# FDR postprocess interrupted @ %s", ctime(&now) );
			break;
		}
		prec += 1;
		n    += 1;
	}
}


/***************************************************************************
  * Row selection iterators
  * Each of these 5 methods takes arguments specific to the
  * source of input pairs and uses file static state for output.
  * All but the first delegate to the _analyze function which may be a
  * filter on immediate-mode output or a caching function for FDR.
  */

// BEGIN:RSI

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
		mtm_resort_rowmap( &_matrix, MTM_RESORT_LEXIGRAPHIC );
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
		mtm_resort_rowmap( &_matrix, MTM_RESORT_LEXIGRAPHIC );
		econd = mtm_fetch_by_name( &_matrix, &(pair.r) );
	} else
		errx( -1, MISSING_MSG, right );

	covan_exec( &pair, &covan );

	_emit( &pair, &covan, g_fp_output );

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

		if( g_sigint_received ) {
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
static int /*AALL*/ _analyze_all_pairs() {

	struct feature_pair fpair;

	struct mtm_row_id *lrid, *rrid;

	assert( ! _matrix.lexigraphic_order /* should be row order */ );

	fpair.l.data = _matrix.data;
	lrid         = _matrix.row_map; // may be NULL

	for(fpair.l.offset = 0; 
		fpair.l.offset < _matrix.rows; 
		fpair.l.offset++ ) {

		fpair.l.name = lrid ? lrid->string : "";
		fpair.l.prop = _matrix.prop[ fpair.l.offset ];

		fpair.r.data = fpair.l.data + _matrix.columns;
		rrid = lrid ? lrid + 1 : NULL;

		for(fpair.r.offset = fpair.l.offset+1; 
			fpair.r.offset < _matrix.rows; 
			fpair.r.offset++ ) {

			fpair.r.name = rrid ? rrid->string : "";
			fpair.r.prop = _matrix.prop[ fpair.r.offset ];

			_analyze( &fpair );

			if( g_sigint_received ) {
				time_t now = time(NULL);
				fprintf( stderr, "# main analysis loop interrupted @ %s", ctime(&now) );
				break;
			}

			fpair.r.data += _matrix.columns;
			if( rrid )  rrid += 1;

		} // inner for

		fpair.l.data += _matrix.columns;
		if( lrid ) lrid += 1;
	}
	return 0;
}

// END:RSI


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
			exename, _VER_MAJ, _VER_MIN, _VER_FIX, _VER_TAG, debug_state,
			exename,
			exename,
			TYPE_PARSER_INFER,
			opt_na_regex,
			DEFAULT_COROUTINE,

			arg_min_cell_count,
			arg_min_sample_count,
			arg_min_mixb_count,

			opt_status_mask,
			opt_p_value,
			opt_format,
			MAGIC_STD,
			MAGIC_TCGA,
			_YN(opt_warnings_are_fatal),
			opt_verbosity,
			MAGIC_SUFFIX,
			MAX_CATEGORY_COUNT,
			MTM_MAX_MISSING_VALUES,
			AUTHOR_EMAIL );
	else
		fprintf( fp, USAGE_ABRIDGED,
			exename, _VER_MAJ, _VER_MIN, _VER_FIX, _VER_TAG, debug_state,
			exename,
			opt_p_value,
		  	AUTHOR_EMAIL );
}


int main( int argc, char *argv[] ) {

	int exit_status = EXIT_SUCCESS;

	const char *i_file = NULL;
	const char *o_file = NULL;
	FILE *fp           = NULL;

	/**
	  * Mandatory initializations.
	  */

	opt_na_regex = mtm_default_NA_regex;
	//memset( _filtered, 0, sizeof(_filtered) );
	memset( &_matrix,  0, sizeof(struct mtm_matrix) );

	/**
	 * This auto(de)selects options related to running this executable
	 * under the TCGA web service.
	 */

	{
		const char *ps
			= strstr( argv[0], MAGIC_SUFFIX );
		opt_running_as_webservice = ( NULL != ps )
			&& strcmp(ps, MAGIC_SUFFIX) == 0; // it's really a suffix
		if( opt_running_as_webservice ) {
			opt_header     = false;
			opt_row_labels = true;
		}
	}

	gsl_set_error_handler( _error_handler );

	/**
	 * Argument checks
	 */

	if( argc < 2 ) { // absolute minimum args: <executable name> <input matrix>
		_print_usage( argv[0], stdout, USAGE_SHORT );
		exit( EXIT_SUCCESS );
	}

	/**
	 * Running as a web service implies some changes to default arguments
	 * ...and NOTHING else. For the purpose of minimizing testing paths,
	 * The webservice executable is identical in every respect to the 
	 * command line too and SHOULD BE MAINTAINED SO.
	 */

	if( opt_running_as_webservice ) {
		// TODO: these need to be (re)defined.
		opt_by_name        = true;
		opt_verbosity      = 0;
		opt_p_value        = 1.0; // implies NO filtering.
	}

	do {

		static const char *CHAR_OPTIONS 
#ifdef HAVE_LUA
			= "s:hrt:N:P:n:x:c:DM:p:f:q:v:?X";
#else
			= "hrt:N:P:n:x:DM:p:f:q:v:?X";
#endif

		static struct option LONG_OPTIONS[] = {
#ifdef HAVE_LUA
			{"script",        required_argument,  0,'s'},
#endif
			{"no-header",     no_argument,        0,'h'},
			{"no-row-labels", no_argument,        0,'r'},
			{"type-parser",   required_argument,  0,'t'},
			{"na-regex",      required_argument,  0,'N'},
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
			if( strcmp( MAGIC_STD, optarg ) == 0 )
				_emit = format_standard;
			else
			if( strcmp( MAGIC_TCGA, optarg ) == 0 )
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

	static const char*NAME_STDIN  = "stdin";
	static const char*NAME_STDOUT = "stdout";

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
		i_file = argv[ optind++ ];
		o_file = argv[ optind++ ];
		break;

	default:
		fprintf( stderr,
			"error: too many (%d) positional arguments supplied. Expect 0, 1, or 2.\n", argc - optind );
		while( optind < argc )
			fprintf( stderr, "\t\"%s\"\n", argv[optind++] );
		exit( EXIT_FAILURE );
	}

	if( opt_pairlist_source != NULL 
			&& strcmp( opt_pairlist_source, i_file ) == 0 ) {
		errx( -1, "error: stdin specified (or implied) for both pair list and the input matrix\n" );
	}

	/**
	  * Emit intentions...
	  */

	if( opt_verbosity > 1 ) {

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
			"defaults for: %s\n"
			"       input: %s\n"
			"      select: %s\n"
			"      output: %s\n"
#ifdef _DEBUG
			"       debug: silent: %s\n"
#endif
			,opt_running_as_webservice ? "web service" : "command line",
			i_file,
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

	if( ! opt_running_as_webservice ) {
		if( SIG_ERR == signal( SIGINT, _interrupt ) ) {
			warn( "failed installing interrupt handler\n"
				"\tCtrl-C will terminated gracelessly\n" );
		}
	}

	/**
	  * Choose and open, if necessary, an output stream.
	  */

	g_fp_output 
		= (strcmp( o_file, NAME_STDOUT ) == 0 )
		? stdout
		: fopen( o_file, "w" );

	if( NULL == g_fp_output ) {
		err( -1, "opening output file \"%s\"", o_file );
	}

	if( covan_init( _matrix.columns ) ) {
		err( -1, "error: covan_init(%d)\n", _matrix.columns );
	}

	if( opt_verbosity > 0 ) 
		fprintf( g_fp_output, "# %d rows X %d (data) columns\n", _matrix.rows, _matrix.columns );


	/**
	  * Here the main decision is made regrding feature selection. 
	  * In order of precedence:
	  * 1. single pair
	  * 2. explicit pairs (by name or by offset)
	  * 3. Lua-generated offsets
	  * 4. all-pairs
	  */

	if( opt_single_pair ) {

		_analyze_single_pair( opt_single_pair, opt_row_labels );

	} else
	if( opt_pairlist_source ) {

		FILE *fp
			= strcmp(opt_pairlist_source,NAME_STDIN)
			? fopen( opt_pairlist_source, "r" )
			: stdin;
		
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

		if( USE_FDR_CONTROL ) {
			g_fp_cache = tmpfile();
			if( NULL == g_fp_cache )
				err( -1, "creating a temporary file" );
		}

#ifdef HAVE_LUA
		if( opt_coroutine )
			_analyze_generated_pair_list( _L );
		else
#endif
			_analyze_all_pairs();

		// Post process results if FDR is in effect and the 
		// 1st pass was allowed to complete.
		// (Post-processing involves a repetition of all of the
		//  above loop FOR A SUBSET of the original input.)

		if( USE_FDR_CONTROL ) {
			if( ! g_sigint_received ) {
				fdr_postprocess( g_fp_cache, arg_q_value );
			}
			if( g_fp_cache ) 
				fclose( g_fp_cache );
		}
	}

	if( opt_verbosity > 0 ) {
		//int i;
		fprintf( g_fp_output, "# Filter counts follow:\n" );
#ifdef HAVE_ANALYSIS
		for(i = 0; i < (int)CovarTypeCount; i++ ) {
			// TODO: fprintf( g_fp_output, "# %s %d\n", COVAR_TYPE_STR[i], _filtered[i] );
		}
#endif
	}

	fclose( g_fp_output );

	return exit_status;
}

