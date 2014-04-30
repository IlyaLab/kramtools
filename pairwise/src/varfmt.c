
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/types.h>
#include <regex.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>

#include "stattest.h"

struct feature_pair;

#include "analysis.h"
#include "varfmt.h"
#include "mtmatrix.h"
#include "featpair.h"

/**
  * Following two variables are initialized by emit_config
  * and used at runtime (by emit_exec) to actually produce output.
  */
static int _columns = 0;
static EMITTER_FXN _emitter[ MAX_OUTPUT_COLUMNS ];
static const char *NOTIMPL = "unimplemented";

/**
  * This holds a copy of the format specification provided to emit_build.
  * The pointers in the printf_format array point *into* this buffer.
  * A subset of the actual _emit{Tab|JSON}* methods use these format
  * strings.
  */
static       char *_emitter_specifier = NULL;

static const char *printf_format[6] = {
	"%.3e", // cov
	"%.3e", // cov
	"%.3e", // cov
	"%.3e", // cov
	"%.3e", // uni
	"%.3e"  // uni
};
#define NUM_PRINTF_FORMATS (sizeof(printf_format)/sizeof(const char*))

static void _emitJSONCovCount( EMITTER_SIG ) {}
static void _emitJSONCovStatName( EMITTER_SIG ) {}
static void _emitJSONCovErrBits( EMITTER_SIG ) {}
static void _emitJSONCovSign( EMITTER_SIG ) {}
static void _emitJSONCovStatValue( EMITTER_SIG ) {}
static void _emitJSONCovProb( EMITTER_SIG ) {}
static void _emitJSONCovLogProb( EMITTER_SIG ) {}
static void _emitJSONCovExtra( EMITTER_SIG ) {}
static void _emitJSONUniNameL( EMITTER_SIG ) {}
static void _emitJSONUniNameR( EMITTER_SIG ) {}
static void _emitJSONUniOffsetL( EMITTER_SIG ) {}
static void _emitJSONUniOffsetR( EMITTER_SIG ) {}
static void _emitJSONUniFeatureClassL( EMITTER_SIG ) {}
static void _emitJSONUniFeatureClassR( EMITTER_SIG ) {}
static void _emitJSONUniPreprocL( EMITTER_SIG ) {}
static void _emitJSONUniPreprocR( EMITTER_SIG ) {}
static void _emitJSONUniUnusedL( EMITTER_SIG ) {}
static void _emitJSONUniUnusedR( EMITTER_SIG ) {}
static void _emitJSONUniStatNameL( EMITTER_SIG ) {}
static void _emitJSONUniStatNameR( EMITTER_SIG ) {}
static void _emitJSONUniStatValueL( EMITTER_SIG ) {}
static void _emitJSONUniStatValueR( EMITTER_SIG ) {}
static void _emitJSONUniProbL( EMITTER_SIG ) {}
static void _emitJSONUniProbR( EMITTER_SIG ) {}
static void _emitJSONUniExtraL( EMITTER_SIG ) {}
static void _emitJSONUniExtraR( EMITTER_SIG ) {}

static void _emitTabCovCount( EMITTER_SIG ) {
	fprintf( fp, "%d", covan->result.sample_count );
}
static void _emitTabCovStatName( EMITTER_SIG ) {
	fputs( covan->result.name, fp );
}
static void _emitTabCovErrBits( EMITTER_SIG ) {
	fprintf( fp, "%02X", covan->status );
}
static void _emitTabCovSign( EMITTER_SIG ) {
	fprintf( fp, printf_format[0], covan->sign );
}
static void _emitTabCovStatValue( EMITTER_SIG ) {
	fprintf( fp, printf_format[1], covan->result.value );
}
static void _emitTabCovProb( EMITTER_SIG ) {
	fprintf( fp, printf_format[2], covan->result.probability );
}
static void _emitTabCovLogProb( EMITTER_SIG ) {
	const double d
		= -log10(covan->result.probability);
	fprintf( fp, printf_format[3], d );
}
static void _emitTabCovExtra( EMITTER_SIG ) {
	fputs( NOTIMPL, fp );
}


static void _emitTabUniNameL( EMITTER_SIG ) {
	fputs( pair->l.name ? pair->l.name : "?", fp );
}
static void _emitTabUniNameR( EMITTER_SIG ) {
	fputs( pair->r.name ? pair->r.name : "?", fp );
}
static void _emitTabUniOffsetL( EMITTER_SIG ) {
	fprintf( fp, "%d",  pair->l.offset );
}
static void _emitTabUniOffsetR( EMITTER_SIG ) {
	fprintf( fp, "%d",  pair->r.offset );
}

static void _emitTabUniFeatureClassL( EMITTER_SIG ) {
#ifndef _UNIT_TEST_VARFMT
	fputs( mtm_sclass_name( covan->stat_class.left ), fp );
#endif
}
static void _emitTabUniFeatureClassR( EMITTER_SIG ) {
#ifndef _UNIT_TEST_VARFMT
	fputs( mtm_sclass_name( covan->stat_class.right ), fp );
#endif
}
static void _emitTabUniPreprocL( EMITTER_SIG ) {
	fputs( NOTIMPL, fp );
}
static void _emitTabUniPreprocR( EMITTER_SIG ) {
	fputs( NOTIMPL, fp );
}
static void _emitTabUniUnusedL( EMITTER_SIG ) {
	fprintf( fp, "%d", covan->waste[0].unused );
}
static void _emitTabUniUnusedR( EMITTER_SIG ) {
	fprintf( fp, "%d", covan->waste[1].unused );
}
static void _emitTabUniStatNameL( EMITTER_SIG ) {
	fputs( covan->waste[0].result.name, fp );
}
static void _emitTabUniStatNameR( EMITTER_SIG ) {
	fputs( covan->waste[1].result.name, fp );
}
static void _emitTabUniStatValueL( EMITTER_SIG ) {
	fprintf( fp, printf_format[4], covan->waste[0].result.value );
}
static void _emitTabUniStatValueR( EMITTER_SIG ) {
	fprintf( fp, printf_format[4], covan->waste[1].result.value );
}
static void _emitTabUniProbL( EMITTER_SIG ) {
	fprintf( fp, printf_format[5], covan->waste[0].result.probability );
}
static void _emitTabUniProbR( EMITTER_SIG ) {
	fprintf( fp, printf_format[5], covan->waste[1].result.probability );
}
static void _emitTabUniExtraL( EMITTER_SIG ) {
	fputs( NOTIMPL, fp );
}
static void _emitTabUniExtraR( EMITTER_SIG ) {
	fputs( NOTIMPL, fp );
}

/**
  * Notice I'm using a restricted printf formatting string.
  * Pad spaces are not allowed, mostly to simplify format parsing.
  */
#define ALLOWED_PRINTF_FORMAT "(%[-#0+]?([1-9][0-9]*)?(.[0-9]*)?[eEfFgG])?"


/**
  * This struct is the "inventory" of all possible output emitters.
  * The _emitter array is actually used at runtime. During program
  * initialization function pointers are copied under the direction
  * of the format specifier from the following struct array into
  * the _emitter array.
  */
struct EmitterDescriptor {
	const char *name;
	const char *pattern;
	regex_t     re;
	EMITTER_FXN emitter[4];
	int         microformat;
} inventory[] = {
	{
		name: "sample count",
		pattern:"^c(ount)?$",
		emitter:{
			_emitTabCovCount,
			NULL,
			_emitJSONCovCount,
			NULL},
		microformat:-1
	},{
		name:"statistic name",
		pattern:"^st(at)?$",
		emitter:{
			_emitTabCovStatName,
			NULL,
			_emitJSONCovStatName,
			NULL},
		microformat:-1
	},{
		name:"test error",
		pattern:"^e(rror)?$",
		emitter:{
			_emitTabCovErrBits,
			NULL,
			_emitJSONCovErrBits,
			NULL},
		microformat:-1
	},{
		name:"(correlation) sign",
		pattern:"^si(gn)?" ALLOWED_PRINTF_FORMAT "$",
		emitter:{
			_emitTabCovSign,
			NULL,
			_emitJSONCovSign,
			NULL},
		microformat:0
	},{
		name:"statistic value",
		pattern:"^v(alue)?" ALLOWED_PRINTF_FORMAT "$",
		emitter:{
			_emitTabCovStatValue,
			NULL,
			_emitJSONCovStatValue,
			NULL},
		microformat:1
	},{
		name:"p-value",
		pattern:"^p(rob)?" ALLOWED_PRINTF_FORMAT "$",
		emitter:{
			_emitTabCovProb,
			NULL,
			_emitJSONCovProb,
			NULL},
		microformat:2
	},{
		name:"-log10(p-value)",
		pattern:"^P(rob)?" ALLOWED_PRINTF_FORMAT "$",
		emitter:{
			_emitTabCovLogProb,
			NULL,
			_emitJSONCovLogProb,
			NULL},
		microformat:3
	},{
		name:"extra",
		pattern:"^e?x(tra)?$",
		emitter:{
			_emitTabCovExtra,
			NULL,
			_emitJSONCovExtra,
			NULL},
		microformat:-1
	},{
		name:"feature",
		pattern:"^[<>]{1,2}f(eature)?$",
		emitter:{
			_emitTabUniNameL,
			_emitTabUniNameR,
			_emitJSONUniNameL,
			_emitJSONUniNameR},
		microformat:-1
	},{
		name:"offset (0-based matrix row)",
		pattern:"^[<>]{1,2}o(ffset)?$",
		emitter:{
			_emitTabUniOffsetL,
			_emitTabUniOffsetR,
			_emitJSONUniOffsetL,
			_emitJSONUniOffsetR},
		microformat:-1
	},{
		name:"statistical class",
		pattern:"^[<>]{1,2}cl(ass)?$",
		emitter:{
			_emitTabUniFeatureClassL,
			_emitTabUniFeatureClassR,
			_emitJSONUniFeatureClassL,
			_emitJSONUniFeatureClassR},
		microformat:-1
	},{
		name:"preprocessing",
		pattern:"^[<>]{1,2}pre(proc)?$",
		emitter:{
			_emitTabUniPreprocL,
			_emitTabUniPreprocR,
			_emitJSONUniPreprocL,
			_emitJSONUniPreprocR},
		microformat:-1
	},{
		name:"unused sample count",
		pattern:"^[<>]{1,2}u(nused)?$",
		emitter:{
			_emitTabUniUnusedL,
			_emitTabUniUnusedR,
			_emitJSONUniUnusedL,
			_emitJSONUniUnusedR},
		microformat:-1
	},{
		name:"statistic name",
		pattern:"^[<>]{1,2}s(tat)?$",
		emitter:{
			_emitTabUniStatNameL,
			_emitTabUniStatNameR,
			_emitJSONUniStatNameL,
			_emitJSONUniStatNameR},
		microformat:-1
	},{
		name:"statistic value",
		pattern:"^[<>]{1,2}v(alue)?$",
		emitter:{
			_emitTabUniStatValueL,
			_emitTabUniStatValueR,
			_emitJSONUniStatValueL,
			_emitJSONUniStatValueR},
		microformat:4
	},{
		name:"p-value",
		pattern:"^[<>]{1,2}pro(bability)?" ALLOWED_PRINTF_FORMAT "$",
		emitter:{
			
			_emitTabUniProbL,
			_emitTabUniProbR,
			_emitJSONUniProbL,
			_emitJSONUniProbR},
		microformat:5
	},{
		name:"extra",
		pattern:"^[<>]{1,2}e?x(tra)?$",
		emitter:{
			_emitTabUniExtraL,
			_emitTabUniExtraR,
			_emitJSONUniExtraL,
			_emitJSONUniExtraR},
		microformat:-1
	},{
		name:NULL,
		pattern:NULL,
		emitter:{NULL,NULL,NULL,NULL},
		microformat:-1
	},
};
#define NUM_COV_EMITTERS (8)
#define NUM_UNI_EMITTERS (9)

static void _fini() {

	struct EmitterDescriptor *d = inventory;
	while( d->name ) {
		regfree( &(d->re) );
		d++;
	}
	if( _emitter_specifier )
		free( _emitter_specifier );
}


/**
  * Called implicitly and lazily within emit_config.
  * This function does not rely on any user-supplied argument
  * so should never fail (post debugging).
  */
static int _init() {

	struct EmitterDescriptor *d = inventory;

	assert( NUM_COV_EMITTERS + NUM_UNI_EMITTERS + 1 
		== sizeof(inventory)/sizeof(struct EmitterDescriptor) );

	while( d->name ) {
#ifdef _UNIT_TEST_VARFMT
		int econd =
#endif
			regcomp( &(d->re), d->pattern, REG_EXTENDED | REG_NOSUB | REG_NEWLINE );
#ifdef _UNIT_TEST_VARFMT
		if( econd ) {
			const size_t required 
				= regerror( econd, &(d->re), NULL, 0 );
			char *buf 
				= (char *)alloca( required );
			regerror( econd, &(d->re), buf, required );
			fprintf( stderr, "%ld: %s\n", d - inventory, buf );
			break;
		}
#endif
		d++;
	}
	atexit( _fini );
	return 0;
}


/**
  * If a specifier is encountered that cannot be matched with 
  * valid patterns, this specifier is returned. Otherwise, NULL.
  */
const char *emit_config( const char *specifier_sequence, int format ) {

	struct EmitterDescriptor *d;
	char *specifier, *pc;
	const int BASE_FORMAT = 2*format;

	pc = _emitter_specifier = strdup( specifier_sequence );

	assert( 0 <= format && format < 2 );

	if( format == FORMAT_JSON ) {
		return "JSON implementation incomplete";
	}

	if( _init() ) 
		return "bug"; // should NEVER happen in Release build.

	while( *pc && _columns < MAX_OUTPUT_COLUMNS ) {

		// Find start of specifier; skip leading space.

		while( *pc &&   isspace(*pc) ) ++pc;

		if( *pc == 0 ) break; // eof of string found; we're done.

		specifier = pc;

		// Find end of specifier.

		while( *pc && ! isspace(*pc) ) ++pc;

		// Terminate and move pointer forward IFF we have NOT
		// encountered EOS.

		if( *pc )
			*pc++ = '\0';

		// Search for a matching regex.

		d = inventory;
		while( d->name ) {
			if( regexec( &(d->re), specifier, 0, NULL, 0 ) == 0 ) {

				if( d->emitter[1] == NULL ) {

					// EmitterDescriptors with NULL in emitter[1] (and 
					// emitter[3]) are covariate emitters and don't need
					// left/right differentiation...

					_emitter[ _columns++ ] = d->emitter[ BASE_FORMAT + 0 ];

				} else { 
		
					// Differentiate between left/right features.

					const char *pc 
						= specifier; // Yes, intentionally mask other pc!
					while( ( *pc == '<' || *pc == '>' )
						&& (_columns < MAX_OUTPUT_COLUMNS ) ) {

						const int LR = *pc == '<' ? 0 : 1;
						_emitter[ _columns++ ] 
							= d->emitter[ BASE_FORMAT + LR ];
						++ pc;
					}
				}
				// Some emitters support an optional microformat.
				if( d->microformat >= 0 ) {
					const char *fmt = strchr( specifier, '%' );
					if( fmt )
						printf_format[ d->microformat ] = fmt;
				}
#ifdef _UNIT_TEST_VARFMT
				fprintf( stdout, "%s\n", d->name );
#endif
				break;
			}
			d++;
		}
		if( d->name == NULL ) {
			break;
		}
	}

	// If we've ever reached the sentinel at the end of the inventory
	// we failed to match an emitter specifier...error!

	if( d->name == NULL )
		return specifier;

	return NULL;
}


void emit_exec( EMITTER_SIG ) {
	int i;
	char *sep = "";
	for(i = 0; i < _columns; i++ ) {
		fputs( sep, fp );
		_emitter[i]( pair, covan, fp );
		sep = "\t";
	}
	fputc( '\n', fp );
}


#ifdef _UNIT_TEST_VARFMT
int main( int argc, char *argv[] ) {
	if( argc > 1 ) {
		emit_config( argv[1], argc > 2 ? atoi(argv[2]) : 0 );
	} else
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
#endif

