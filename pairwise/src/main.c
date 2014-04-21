
/**
  * This file serves 3 primary functions:
  * 1) finding and prepping the input matrix for processing.
  * 2) iterating through the input matrix' rows either by name or 
  *    using a possibly complex series of range specifications
  *    and for each pair calling a high-level analysis method
  * 3) routing the results to either
  *    a) immediate output using a specified formatter or
  *    b) a cache to be post-processed (for FDR control) and subsequently
  *       output.
  *
  *     
  *
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

//#include "format.h"
#include "mtmatrix.h"
#include "stattest.h"
#include "analysis.h"

#ifdef HAVE_LUA
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#endif

/***************************************************************************
 * externs
 */

extern int get_base10_ints( FILE *fp, int *index, int n );
extern int sclass_by_prefix( const char *token );

/***************************************************************************
 * Globals & statics
 */

FILE *g_fp_output = NULL;
FILE *g_fp_cache = NULL;

size_t g_COLUMNS = 0;
/*
void (*g_format)( 
		unsigned lf, 
		unsigned rf,
		unsigned status ) = format_tcga;
*/

#define GLOBAL
GLOBAL unsigned  arg_min_cell_count   = 5;
GLOBAL unsigned  arg_min_mixb_count   = 1;
GLOBAL unsigned  arg_min_sample_count = 2; // < 2 NEVER makes sense
GLOBAL bool g_sigint_received          = false; // fdr.d needs access
#undef  GLOBAL

static const char *TYPE_PARSER_INFER   = "auto";
static double arg_q_value              = 0.0;
#define USE_FDR_CONTROL (arg_q_value > 0.0)

static const char *opt_script          = NULL;
static bool        opt_header          = true;
static bool        opt_row_labels      = true;
static const char *opt_type_parser     = NULL;
static const char *opt_na_regex        = "[Nn][Aa][Nn]?";

static const char *opt_single_pair     = NULL;

static const char *opt_pairlist_source = NULL;
static bool        opt_by_name         = false;
static const char *DEFAULT_COROUTINE   = "pair_generator";
static const char *opt_coroutine       = NULL;
static bool        opt_dry_run         = false;

static double      opt_p_value         = 1.0;
static const char *opt_format          = "foo";

static bool     opt_warnings_are_fatal = false;

static int         arg_verbosity       = 0;

static bool  arg_webservice       = false;

#ifdef _DEBUG
/**
  * This debug-only option is for the purpose of full-pass testing
  * without loads of irrelevant output.
  */
static bool  dbg_exhaustive = false;

/**
  * Emit all results INCLUDING those with degeneracy.
  */

static bool  dbg_silent     = false;
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
static struct mtmatrix _matrix;

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
 * B) NC2 or Lua generated pair lists
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
#define ANALYSIS_FN_SIG const struct mt_row_pair *pair

typedef void (*ANALYSIS_FN)( ANALYSIS_FN_SIG );



/***************************************************************************
 * Encapsulates all the decision making regarding actual emission of
 * results.
 */
static void _filter( ANALYSIS_FN_SIG ) {

	struct CovariateAnalysis covan;
	memset( &covan, 0, sizeof(covan) );
	covan_exec( pair, &covan );

#ifdef _DEBUG
	if( ! dbg_silent ) {
#endif

		if(
#ifdef _DEBUG
			dbg_exhaustive || 
#endif
			(( covan.result.probability <= opt_p_value ) && ( covan.status == 0 ) ) 

			) {
			//g_format( a, b, status );
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
  * Note that this method does NOT take a cache FILE* like fdr_postprocess
  * ONLY because it must follow a prototype that makes it interchangeable
  * with the _emit static in main. 
  * TODO: Revisit this.
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
		//g_format( prec->a, prec->b, status );

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
  * Each of these (currently 4) methods takes arguments specific to the
  * source of input pairs and uses file static state for output.
  */

static int /*ANAM*/ _analyze_named_pair_list( FILE *fp ) {

	struct mt_row_pair fpair;

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

		if( mtm_fetch_by_name( &_matrix, &fpair ) != 2 ) {

			fprintf( stderr,
				"error: one or both of...\n"
				"\t1) %s\n"
				"\t2) %s\n"
				"\t...not found.\n",
				fpair.l.name,
				fpair.r.name );
			if( opt_warnings_are_fatal ) 
				break;
			else
				continue; // no reason we -can't- continue
		}

		_analyze( &fpair );
	}

	if( left )
		free( left );

	return 0;
}


static int /*ANUM*/ _analyze_pair_list( FILE *fp ) {

	struct mt_row_pair fpair;

	int arr[2];

	while( get_base10_ints( fp, arr, 2 ) == 2 ) {

		fpair.l.offset  = arr[0];
		fpair.r.offset = arr[1];

		if( mtm_fetch_by_offset( &_matrix, &fpair ) != 2 ) {
			fprintf( stderr, 
				"error: one of row indices (%d,%d) not in [0,%d)\n"
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
		}
		
		_analyze( &fpair );
	}
	return 0;
}


#ifdef HAVE_LUA
static int /*ALUA*/ _analyze_generated_pair_list( lua_State *state ) {

	struct mt_row_pair fpair;

	while( false ) {
		_analyze( &fpair );
		if( g_sigint_received ) {
			time_t now = time(NULL);
			fprintf( stderr, "# main analysis loop interrupted @ %s", ctime(&now) );
			break;
		}
	}

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

	struct mt_row_pair fpair;

	struct row_id *lrid, *rrid;

	assert( ! _matrix.lexigraphic_order /* should be row order */ );

	fpair.l.data = _matrix.data;
	lrid         = _matrix.row_map;

	for(fpair.l.offset = 0; 
		fpair.l.offset < _matrix.rows; 
		fpair.l.offset++ ) {

		fpair.l.name = lrid ? lrid->string : "";
		fpair.l.prop = _matrix.prop[ fpair.l.offset ];

		fpair.r.data = fpair.l.data + _matrix.columns;
		rrid = lrid + 1;

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


/**
  * <source> may be literal Lua code or a filename reference.
  */
static int _load_script( const char *source, lua_State *state ) {

	struct stat info;
	return (access( source, R_OK ) == 0) && (stat( source, &info ) == 0)
		? luaL_dofile(   state, source )
		: luaL_dostring( state, source );
}


/***************************************************************************
  * Online help
  */


static const char *_YN( bool y ) {
	return y ? "yes" : "no";
}

static void _print_usage( const char *exename, FILE *fp, bool exhaustive ) {

	extern const char *USAGE_UNABRIDGED;
	//extern const char *USAGE_ABBREV;

	if( exhaustive || ! exhaustive )
		fprintf( fp, USAGE_UNABRIDGED,
			exename, _VER_MAJ, _VER_MIN, _VER_FIX,
#ifdef _DEBUG
			"(DEBUG)",
#else
			"",
#endif
			exename,
			TYPE_PARSER_INFER,
			opt_na_regex,
			DEFAULT_COROUTINE,

			arg_min_cell_count,
			arg_min_sample_count,
			arg_min_mixb_count,

			opt_p_value,
			opt_format,
			_YN(opt_warnings_are_fatal),
			arg_verbosity,
			_MAX_CATEGORIES );
}


int main( int argc, char *argv[] ) {

	int exit_status = EXIT_SUCCESS;

	const char *i_file = NULL;
	const char *o_file = NULL;
	FILE *fp           = NULL;

#ifdef HAVE_LUA
	lua_State *L = NULL;
#endif

	/**
	  * Mandatory initializations.
	  */

	//memset( _filtered, 0, sizeof(_filtered) );
	memset( &_matrix,  0, sizeof(struct mtmatrix) );

	/**
	 * This auto(de)selects options related to running this executable
	 * under the TCGA web service.
	 */

	{
		const char *ps
			= strstr( argv[0], "-www" );
		arg_webservice 
			= ( NULL != ps ) && strlen(ps) == 4; // confirm it's at the end.
	}

	gsl_set_error_handler( _error_handler );

	/**
	 * Argument checks
	 */

	if( argc < 2 ) { // absolute minimum args: <executable name> <input matrix>
		_print_usage( argv[0], stdout, false );
		exit(0);
	}

	/**
	 * Running as a web service implies some changes to default arguments
	 * ...and NOTHING else. For the purpose of minimizing testing paths,
	 * The webservice executable is identical in every respect to the 
	 * command line too and SHOULD BE MAINTAINED SO.
	 */

	if( arg_webservice ) {
		opt_by_name        = true;
		arg_verbosity      = 0;
		opt_p_value        = 1.0; // implies NO filtering.
	}

	do {

		static const char *CHAR_OPTIONS 
			= "s:hrt:N:P:n:x:c:DM:p:f:q:v:?X";

		static struct option LONG_OPTIONS[] = {
			{"script",        required_argument,  0,'s'},
			{"no-header",     no_argument,        0,'h'},
			{"no-row-labels", no_argument,        0,'r'},
			{"type-parser",   required_argument,  0,'t'},
			{"na-regex",      required_argument,  0,'N'},
			{"pair",          required_argument,  0,'P'},

			{"by-name",       required_argument,  0,'n'},
			{"by-index",      required_argument,  0,'x'},
			{"coroutine",     required_argument,  0,'c'},
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
				fprintf( stderr, 
					"Seriously...%d samples is acceptable?\n"
					"I don't think so...\n", 
					arg_min_sample_count );
				exit(-1);
			}
			break;
		////////////////////////////////////////////////////////////////////

		case 'p':
			opt_p_value = atof( optarg );
			if( ! ( 0 < opt_p_value ) ) {
				fprintf( stderr,
					"error: specified p-value %.3f will preclude all output.\n",
					opt_p_value );
				abort();
			} else
			if( ! ( opt_p_value < 1.0 ) ) {
				fprintf( stderr,
					"warning: p-value %.3f will filter nothing.\n"
					"\tIs this really what you want?\n", 
					opt_p_value );
			}
			break;

		case 'f': // format
			if( strncmp("std",optarg, 3 ) == 0 )
				;// TODO: g_format = format_standard;
			else
			if( strncmp("short",optarg, 5 ) == 0 )
				;// TODO: g_format = format_abbreviated;
			break;

		case 'q':
			arg_q_value = atof( optarg );
			break;

		case 'v': // verbosity
			arg_verbosity = atoi( optarg );
			break;

		case '?': // help
			_print_usage( argv[0], stderr, false );
			exit(0);
			break;

		case 'X': // help
			_print_usage( argv[0], stderr, true );
			exit(0);
			break;

#ifdef _DEBUG
		case 258:
			dbg_exhaustive = strchr( optarg, 'X' ) != NULL;
			dbg_silent     = strchr( optarg, 'S' ) != NULL;
			break;
#endif
		case -1: // ...signals no more options.
			break;
		default:
			printf ("error: unknown option: %c\n", c );
			exit(-1);
		}
		if( -1 == c ) break;

	} while( true );

#ifdef HAVE_LUA

	/**
	  * Load the Lua from file or command line before anything else.
	  */

	if( opt_script ) {

		L = luaL_newstate();
		if( L ) {
			luaL_openlibs(L);
			if( _load_script( opt_script, L ) != LUA_OK ) {
				errx( -1, "failed executing Lua \"%s\": %s", 
					opt_script, lua_tostring(L,-1) );
				lua_close(L);
			}
		} else
			errx( -1, "failed creating Lua statespace" );

		// Additionally, if a Lua script was provided but no coroutine
		// was specified, see whether a function with the default coroutine
		// name exists...

		if( NULL == opt_coroutine ) {

			lua_getglobal( L, DEFAULT_COROUTINE );

			if( ! lua_isnil( L, -1 ) )
				opt_coroutine = DEFAULT_COROUTINE;

			// Non-existence of DEFAULT_COROUTINE is not an error;
			// we assume user does NOT intend Lua to generate the pairs.

			lua_pop( L, 1 );
		}

		/**
		  * Load Lua script and execute a dry-run if applicable.
		  * This is in advance of other argument validation since other
		  * arguments will be ignored...
		  */

		if( opt_dry_run ) {
			int isnum, lua_status = LUA_YIELD;
			const char *coroutine
				= opt_coroutine 
				? opt_coroutine 
				: DEFAULT_COROUTINE;
			lua_getglobal( L, coroutine );
			if( lua_isnil( L, -1 ) )
				errx( -1, "%s not defined (in Lua's global namespace)", 
					coroutine );
			while( lua_status == LUA_YIELD ) {
				lua_status = lua_resume( L, NULL, 0 );
				if( lua_status <= LUA_YIELD /* OK == 0, YIELD == 1*/ ) {
					const int b = lua_tonumberx( L, -1, &isnum );
					const int a = lua_tonumberx( L, -2, &isnum );
					lua_pop( L, 2 );
					fprintf( stdout, "%d %d\n", a, b );
				} else
					fputs( lua_tostring(L,-1), stderr );
			}

			exit(0);
		}
	}

#endif

	/**
	  * Catch simple argument inconsistencies and fail early!
	  */

	if( opt_single_pair ) {
		if( USE_FDR_CONTROL ) {
			fprintf( stderr, "warning: FDR is senseless on a single pair.\n" );
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

	switch( argc - optind ) {

	case 0: // input MUST be stdin, output stdout
		i_file = "stdin";
		o_file = "stdout";
		break;

	case 1:
		if( access( argv[ optind ], R_OK ) == 0 ) {
			i_file = argv[ optind++ ];
			o_file = "stdout";
		} else {
			i_file = "stdin";
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
		abort();
	}

	if( opt_pairlist_source != NULL 
			&& strcmp( opt_pairlist_source, i_file ) == 0 ) {
		fprintf( stderr, 
			"error: stdin specified (or implied) for both pair list and the input matrix\n" );
		abort();
	}

	/**
	  * Load the input matrix.
	  */

	fp = strcmp( i_file, "stdin" )
		? fopen( i_file, "r" )
		: stdin;
	if( fp ) {

		const unsigned int FLAGS
			= ( opt_header ? MATRIX_HAS_HEADER : 0 ) 
			| ( opt_row_labels ? MATRIX_HAS_ROW_NAMES : 0 )
			| 1 /*VERBOSITY_MASK*/;

		const int err
			= mtm_parse( fp,
				FLAGS,
				opt_na_regex,
				_MAX_CATEGORIES,
				sclass_by_prefix,
				&_matrix );
		if( ! err ) 
			atexit( _freeMatrix );
		fclose( fp );
	}

	if( opt_dry_run ) { // a second possible 
		exit(0);
	}

	if( ! arg_webservice ) {
		if( SIG_ERR == signal( SIGINT, _interrupt ) ) {
			fprintf( stderr,
				"warning: failed installing interrupt handler\n"
				"\tCtrl-C will terminated gracelessly\n" );
		}
	}

	g_fp_output 
		= strlen(o_file) > 0 
		? fopen( o_file, "w" ) 
		: stdout;

	if( covan_init( _matrix.columns ) ) {
		fprintf( stderr, "error: covan_init(%d)\n", _matrix.columns );
		exit(-1);
	}

	if( arg_verbosity > 0 ) 
		fprintf( g_fp_output, "# %d rows X %d (data) columns\n", _matrix.rows, _matrix.columns );

	if( opt_single_pair ) {

	} else
	if( opt_pairlist_source ) {

		FILE *fp
			= strcmp(opt_pairlist_source,"stdin")
			? fopen( opt_pairlist_source, "r" )
			: stdin;
		
		if( fp ) {
			int err
				= opt_by_name 
				? _analyze_named_pair_list( fp )
				: _analyze_pair_list( fp );
			fclose( fp );
		} else
			fprintf( stderr, "error" );

	} else {

		if( USE_FDR_CONTROL ) {
			g_fp_cache = tmpfile();
			if( NULL == g_fp_cache ) {
				perror( "creating a temporary file" );
				abort();
			}
		}

#ifdef HAVE_LUA
		if( opt_coroutine )
			_analyze_generated_pair_list( L );
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

	if( arg_verbosity > 0 ) {
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

