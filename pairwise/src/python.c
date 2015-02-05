
#include <Python.h>

#ifdef HAVE_LUA
#error Python and Lua dependencies are mutually exclusive
#endif


extern int mtm_sclass_by_prefix( const char *token );


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



/**
  * All the following is for a very specific ISB-only use case related to a
  * web app developed by Dick Kreisberg et al.
  * ../test/testpyext.py exercises this, but that's your only guide.
  * There is no further documentation and none planned for this.
  */

static PyObject * _run( PyObject *self, PyObject *args ) {

	PyObject * ret = Py_None;
	FILE *fp_input = NULL;
	PyObject *fobj[2];

	if( ! PyArg_ParseTuple( args, "OO", fobj+0, fobj+1 ) ) {
		// PyArg_ParseTuple has already raised appropriate exception.
		return NULL;
	}
	if( ! PyFile_Check( fobj[0] ) ) {
		PyErr_SetString( PyExc_TypeError, "1st argument is not a File object" );
		return NULL;
	}

	if( ! PyFile_Check( fobj[1] ) ) {
		PyErr_SetString( PyExc_TypeError, "2nd argument is not a File object" );
		return NULL;
	}

	 fp_input  = PyFile_AsFile( fobj[0] );
	_fp_output = PyFile_AsFile( fobj[1] );

	if( fp_input != NULL && _fp_output != NULL ) {
	
		int err = 0;

		PyFile_IncUseCount((PyFileObject*)fobj[0]);
		PyFile_IncUseCount((PyFileObject*)fobj[1]);
		Py_BEGIN_ALLOW_THREADS

		/**
		  * Load the input matrix.
		  */

		err	= mtm_parse( fp_input,
				MTM_MATRIX_HAS_ROW_NAMES | MTM_MATRIX_HAS_HEADER,
				mtm_default_NA_regex,
				MAX_CATEGORY_COUNT,
				mtm_sclass_by_prefix,
				&_matrix );

		Py_END_ALLOW_THREADS
		PyFile_DecUseCount((PyFileObject*)fobj[1]);
		PyFile_DecUseCount((PyFileObject*)fobj[0]);

		// Python interpreter is locked (out) for remainder of 
		// pairwise execution...

		if( err == 0 ) {

			if( covan_init( _matrix.columns ) == 0 ) {

				if( _analyze_all_pairs() ) {
					PyErr_SetString( PyExc_KeyboardInterrupt, "interruption" );
					// ...only way _analyze_all_pairs returns non-zero.
					ret = NULL;
				}

				covan_fini();

			} else {
				PyErr_SetString( PyExc_MemoryError,
						"pre-allocating analysis buffers" );
				ret = NULL;
			}

			_matrix.destroy( &_matrix );

		} else {

			static char buf[ 40 ];
			// TODO: Improve error reporting in following.
			sprintf( buf, "failed loading matrix (libmtm error %d)", err );
			PyErr_SetString( PyExc_RuntimeError, buf );
			ret = NULL;
		}

	} else {

		PyErr_SetString( PyExc_RuntimeError, "PyFile_AsFile returned NULL" );
		ret = NULL;
	}

	return ret;
}

static PyMethodDef methods[] = {
	{"run",_run,METH_VARARGS},
	{NULL,NULL},
};

PyMODINIT_FUNC initpairwise(void) {
	PyObject *m
		= Py_InitModule( "pairwise", methods );
	_jit_initialization();
}

#endif
