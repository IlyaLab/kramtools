
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <err.h>

#include "fnv/fnv.h"

#include "strset.h"
#include "mtmatrix.h"
#include "feature.h"
#include "toktype.h"
#include "mterror.h"
#include "mtsclass.h"
#include "specialc.h"

extern int mtm_sclass_by_prefix( const char *token );
extern int cardinality(
		const unsigned int *buf, int len, int largest_of_interest, const unsigned int NA );

const char *mtm_default_NA_regex = "^[Nn][Aa][Nn]?$";
static const char *_BUG = "bug at %s:%d";

/**
  * Using separated fields, column counting is trivial: there will
  * AlWAYS be one more column than there are separators.
  *
  * If the line contains only separators, all fields are empty!
  *
  * Ihis, in particular, allows for the (R-style) of an empty
  * left-justified field in the header line.
  */
int feature_count_fields( const char *pc, const char SEP ) {

	int n = 0;
	while( *pc ) if( *pc++ == SEP ) n++;
	return ++n; // ..."++" for the implicit last token.
}


/**
  * This is atomic: if it doesn not entirely succeed, it deallocates/
  * frees whatever -was- allocated and returns CLEANLY.
  */
int feature_alloc_encode_state( struct feature *f ) {

	assert( sizeof(mtm_int_t) == sizeof(mtm_fp_t) );
	assert( f->length > 1 );

	/**
	  * Initialize the "missing data" detector.
	  */
	if( toktype_init( f->missing_data_regex 
			? f->missing_data_regex : mtm_default_NA_regex ) ) {
		return MTM_E_SYS;
	}

	// Buffers

	f->buf.num = calloc( f->length, sizeof(mtm_fp_t) );
	if( NULL == f->buf.num ) {
		goto failure;
	}

	// Hash table for counting categorical variables' categories.
	// Failure can only be due to out-of-memory condition.

	// Notice: using column count for the string set means szs_insertion
	// will never fail, but it's also potentially wasteful, since it's
	// going to be rounded up to 2^ceil(log_2(column_count)).

	f->category_labels = szs_create( f->length,
			0, // strdup not necessary
			fnv_32_str, FNV1_32_INIT );
	if( f->category_labels == NULL ) {
		goto failure;
	}

	return MTM_OK;

failure:
	feature_free_encode_state( f );

	return MTM_E_IO;
}


/**
  * Convert a string containing some number of fields into binary form,
  * possibly inferring the type of data in the process.
  *
  * Strings are mapped in the order they first appear in the input (so the
  * first string is ALWAYS mapped to 0, etc).
  *
  * ASSUMPTIONS:
  * 1. <line> is NUL-terminated and free of NL and CR characters.
  */
int feature_encode( char * const line,
		struct feature *f,
		struct mtm_descriptor *d ) {

	int stat_class = MTM_STATCLASS_UNKNOWN;
	bool infer_field_type = false;
	int ft, field_type = MTM_FIELD_TYPE_UNK;
	const char *token;
	char *pc = line;

	bool eol = false;
	unsigned int field_count = 0;
	int missing_value_count = 0;

	union {
		mtm_fp_t  f;
		mtm_int_t i;
	} last_value_read;

	memset( d, 0, sizeof(struct mtm_descriptor) );
	assert( szs_count( f->category_labels ) == 0 );

	if( f->expect_row_labels ) {

		// Scan past the row identifier.

		while( *pc && *pc != CHAR_FIELD_SEP ) pc++;

		// ...and make sure there was more than just a row id on the line!

		if( CHAR_FIELD_SEP != *pc ) {
			// A non-empty line with exactly one field is bad format.
			return MTM_E_FORMAT_MATRIX;
		}

		f->label_length = pc - line;
		*pc++ = 0;

		if( f->interpret_prefix ) {
			stat_class
				= f->interpret_prefix( line );
			field_type
				= field_type_from_stat_class( stat_class );
			infer_field_type
				= (__builtin_popcount( field_type ) != 1 );
				// ...either MTM_FIELD_TYPE_UNK or multiple bits.
		}
	}

	while( field_count < f->length && ! eol ) {

		char *endpt = ""; // ESSENTIAL initialization.
		token = pc;

		while( *pc && *pc != CHAR_FIELD_SEP )
			pc++;

		if( token == pc /* we didn't move */ ) {

			++missing_value_count;
			f->buf.cat[ field_count++ ] = NAN_AS_UINT;

			if( *pc++ )
				continue;
			else
				break; // ...because we're already ON eol.
		}

		// We have a non-empty token.
		// Before NUL-terminating see whether it's also the end-of-line.

		if( *pc == 0 )
			eol = true; // We HAVE last field now; preclude loop reentry!
		else
			*pc++ = 0;

		// The "missing data" marker is common to all rows of the matrix
		// so its detection can precede (and preclude) type-dependent ops.

		if( toktype_is_na_marker(token) ) {
			++missing_value_count;
			f->buf.cat[ field_count++ ] = NAN_AS_UINT;
			continue;
		}

		// For every line, we will -not- know the field_type at this point
		// for at most one non-missing token. In other words, almost every
		// time we reach this point field_type is known. Use of the goto thus
		// obviates a constitutive conditional...justifying an ugly goto.
retry:
		switch( field_type ) {

		case MTM_FIELD_TYPE_FLT:
			f->buf.num[ field_count++ ]
				= last_value_read.f
				= strtof( token, &endpt );
			break;

		case MTM_FIELD_TYPE_INT:
			f->buf.cat[ field_count++ ]
				= last_value_read.i
				= strtol( token, &endpt, 0 );
			if( ! ( last_value_read.i < NAN_AS_UINT ) ) {
				return MTM_E_LIMITS;
			}
			break;

		case MTM_FIELD_TYPE_STR:
			if( szs_insert( f->category_labels, token, &(last_value_read.i) ) < 0 )
				errx( -1, _BUG, __FILE__, __LINE__ );
			// Because I'm sizing the hashtable to accomodate a cardinality
			// equal to the sample count, szs_insert cannot fail, or rather
			// if it fails something else is badly broken.
			f->buf.cat[ field_count++ ] = last_value_read.i;
			break;

		default: // ...because field_type is either _UNKNOWN or set valued.
			// This insures that boolean data represented by {0,1} is parsed
			// as integral data, not strings, and thus has its implicit
			// ordering preserved. This was both original code and a bugfix. 
			ft = toktype_infer_narrowest_type( token, NULL );
			assert( __builtin_popcount(ft)==1 /* else infinite loop! */ );
			// If the type was constrained at all (not simply "unknown"),
			// then the inferred type MUST be one of the allowed types...
			if( field_type!=MTM_FIELD_TYPE_UNK && (field_type & ft)==0 ) {
				return MTM_E_FORMAT_FIELD; // It's not an allowed type!
			} else
				field_type = ft;

			goto retry; // yes, goto! See comments above..
		}

		assert( __builtin_popcount( field_type ) == 1 );

		/**
		  * REVISE INFERENCE if necessary:
		  *
		  * If endpt wasn't advanced to the end-of-string the token
		  * does not contain the type it was believed to contain.
		  * This is either because of an overt error or because the stat
		  * class didn't fully constrain token type and the inference from
		  * preceding fields was too narrow. Specifically either:
		  * 1) a floating point line had an integral first value
		  * 2) a string line had a numeric (integral or float) first value
		  * As long as type wasn't completely determined by stat class--that
		  * is, as long as the field_type was inferred--we can revise our
		  * inference in these cases...
		  */

		if( *endpt /* token wasn't entirely consumed */ ) {

			// If the type was dictated rather than inferred we've
			// encountered an error; revision is only possible on
			// inferred types.

			if( infer_field_type ) {

				if( (field_type == MTM_FIELD_TYPE_INT)
					&& (MTM_FIELD_TYPE_FLT == toktype_infer_narrowest_type( token, NULL ) ) ) {
#ifdef _DEBUG
					fputs( "promoting integral line to float\n", stderr );
#endif
					field_type = MTM_FIELD_TYPE_FLT;
					--field_count;
					// Convert all earlier values to float...
					for(int i = 0; i < field_count; i++ ) {
						if( NAN_AS_UINT != f->buf.cat[ i ] ) {
							f->buf.num[i] = (float)f->buf.cat[i]; // in place!
						}
					}
					// ...and reparse current token.
					f->buf.num[ field_count++ ]
						= strtof( token, &endpt );
				} else
				if( (field_type != MTM_FIELD_TYPE_STR)
					&& (MTM_FIELD_TYPE_STR == toktype_infer_narrowest_type( token, NULL ) ) ) {
#ifdef _DEBUG
					fputs( "revising non-string line to string\n", stderr );
#endif
					assert( szs_count( f->category_labels ) == 0 );
					field_type = MTM_FIELD_TYPE_STR;

					// Reparse line up to and including current token.
					// Above inference of eol still holds, and, importantly,
					// we can count on NUL-termination of ALL relevant
					// tokens.

					pc = line;
					field_count = 0;
					do {
						if( szs_insert( f->category_labels, pc, &(last_value_read.i) ) < 0 )
							errx( -1, _BUG, __FILE__, __LINE__ );
						// See comment above re: szs_insert.
						f->buf.cat[ field_count++ ] = last_value_read.i;
						if( pc == token) break;
						pc += strlen(pc)+1;
					} while( true );
					pc += strlen(pc)+1;
					endpt = ""; // to remove the error condition.
				}
			}

			if( *endpt ) // either because the line wasn't a candidate for
				// revision or because it was revised, but the current
				// token -still- wasn't entirely consumed.
				return MTM_E_FORMAT_FIELD;
		}
	}

	if( field_count != f->length ) {
		return MTM_E_FORMAT_MATRIX;
	}

	////////////////////////////////////////////////////////////////////////
	// Now fill out the descriptor
	////////////////////////////////////////////////////////////////////////

	d->missing = missing_value_count;

	// Check for the degeneracy. There is really only one kind of
	// degeneracy, constancy, which can manifest in three different ways:
	// 1. All values are missing, i.e. feature is entirely "NA".
	// 2. Exactly one value is not-missing; obviously one value is constant.
	// 3. Two or more values are non-missing, but they are all equal.
	// The first two are trivial to check for; the 3rd is more expensive...

	if( field_count - missing_value_count < 2 )
		d->constant = 1;

	// ...so we only check for the 3rd if necessary (below).

	switch( field_type ) {

	case MTM_FIELD_TYPE_FLT:
		if( ! d->constant &&
				cardinality( f->buf.cat, field_count, 2, NAN_AS_UINT ) < 2 ) {
			d->constant = 1;
		}
		break;

	case MTM_FIELD_TYPE_STR:
		d->integral    = 1;
		d->categorical = 1;
		d->cardinality = szs_count( f->category_labels );
		if( d->cardinality < 2 )
			d->constant = 1;
		szs_clear( f->category_labels );
		break;

	case MTM_FIELD_TYPE_INT:
		// TODO: Note that ordinal is not properly dealt with yet.
		// If cardinality is low-enough integer data is assumed
		// categorical.
		d->integral = 1;
		if( ! d->constant ) {
			d->cardinality
				= cardinality( f->buf.cat, field_count, f->max_cardinality, NAN_AS_UINT );
			d->categorical
				= d->cardinality <= f->max_cardinality
				? 1
				: 0;
		}
		break;

	case MTM_FIELD_TYPE_UNK:
	default:
		if( missing_value_count < field_count )
			errx( -1, _BUG, __FILE__, __LINE__ );
		// ...but if whole line was empty we may, indeed, not know the type!
	}

	if( MTM_STATCLASS_BOOLEAN == stat_class ) {
		if( d->cardinality > 2 )
			return MTM_E_CARDINALITY;
	}

	return MTM_OK;
}


void feature_free_encode_state( struct feature *f ) {

	// Must be tolerant of -incomplete- initialization!

	if( f->category_labels )
		szs_destroy( f->category_labels );

	if( f->buf.num )
		free( f->buf.num );

	memset( f, 0, sizeof(struct feature) );
}

#ifdef UNIT_TEST_FEATURE

// Areas for testing parsing/encoding:
// 	1. data content
// 	2. matrix structure
// 	3. row labels
// 	4. statistical class prefixes (if present)
// 	5. type inference and statistical class selection in absence
// 		of prefixes
// 
// Most of #1 and #2 and the overall integrity of the filesystem binary
// format can be assessed by round-trip evaluation of the data. That is,
// if a binary matrix is echoed back as text it should be identical to 
// the original input.
// 
// This ideal is precluded since 
// 	1. The string information on categorical data is not preserved in the
// 		encoding of the runtime form.
// 	2. Variability in representation of missing data cannot be reproduced
// 		in the echoed representation.
// It is also complicated by float-point rounding.
// However, if test code .
// 
// Overall structure
// // Line parsing and interpretation is independently testable since
// 	EVERYTHING is scoped to a line. This includes:
// 	stat class and type inference handling
// 	missing data handling
// 	content
// 
// Type inference/statistical class
// 
// 	boolean, categorical, numeric
// 	X
// 	random missing data + missing data representation
// 	X
// 	with(out) row label
// 		X
// 		with(out) statclass prefix
// 
// 	Verify:
// 		Descriptor has correct
// 			flags:
// 				constant
// 				integral
// 				categorical
// 			cardinality
// 			missing count
// 
// 	Forced error conditions to inject:
// 		feature tagged boolean with > 2 levels
// 	Confirm
// 		feature tagged boolean with integer representation preserves rep.

#include <getopt.h>

/**
  * Using separated fields, column counting is trivial: there will
  * AlWAYS be one more column than there are separators.
  *
  * If the line contains only separators, all fields are empty!
  *
  * Ihis, in particular, allows for the (R-style) of an empty
  * left-justified field in the header line.
  */

/**
  * The output defined by this function is intended to be straightforward
  * for the peer Python script to validate in combination with the input
  * options.
  */
static void _emit( const char *line, struct feature *f,
		struct mtm_descriptor *d, FILE *fp ) {

#if 1
	fprintf( fp,
		"DEG\t%s\n"
		"INT\t%s\n"
		"CAT\t%s\n", 
		d->constant    ? "T" : "F",
		d->integral    ? "T" : "F",
		d->categorical ? "T" : "F" );
#else
	if( d->constant    )
		fprintf( fp, "constant\n" );
	if( d->integral    )
		fprintf( fp, "integral\n" );
	if( d->categorical )
		fprintf( fp, "categorical\n" );
#endif

	fprintf( fp, 
		"CAR\t%d\n"
		"N/A\t%d\n"
		"LAB\t%s\n",
		d->cardinality,
		d->missing,
		f->label_length > 0 ? line : "-" );

	if( d->integral ) {
		for(int i = 0; i < f->length; i++ ) {
			fprintf( fp, "%d\n", f->buf.cat[i] );
		}
	} else {
		for(int i = 0; i < f->length; i++ ) {
			fprintf( fp, "%f\n", f->buf.num[i] );
		}
	}
}


int main( int argc, char *argv[] ) {

	int exit_status = EXIT_FAILURE;
	const int FIELD_SEP
		= getenv("SEPARATOR")
		? getenv("SEPARATOR")[0]
		: '\t';

	struct feature f = {
		.length = 0,
		//.buf.num: NULL,
		.expect_row_labels =  true,
		.missing_data_regex = mtm_default_NA_regex,
		.interpret_prefix =   mtm_sclass_by_prefix,
		.max_cardinality =    32,
		.category_labels =    NULL
	};

	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;

	do {
		static const char *CHAR_OPTIONS 
			= "rm:k:v:?";
		static struct option LONG_OPTIONS[] = {

			{"nolabels",   0,0,'r'},
			{"missing",    1,0,'m'},
			{"maxcats",    1,0,'k'},

			{"verbosity",  1,0,'v'},
			{ NULL,        0,0, 0 }
		};

		int opt_offset = 0;
		const int c = getopt_long( argc, argv, CHAR_OPTIONS, LONG_OPTIONS, &opt_offset );
		switch (c) {

		case 'r': f.expect_row_labels  = false;        break;
		case 'm': f.missing_data_regex = optarg;       break;
		case 'k': f.max_cardinality    = atoi(optarg); break;
		case 'v':
			//opt_verbosity = atoi( optarg );
			break;
		case -1: // ...signals no more options.
			break;
		default:
			printf ("error: unknown option: %c\n", c );
			exit(-1);
		}
		if( -1 == c ) break;

	} while( true );

	llen = getline( &line, &blen, stdin );

	if( llen > 0 && line != NULL ) {

		int ec = MTM_OK;
		struct mtm_descriptor d;

		f.length 
			= feature_count_fields( line, FIELD_SEP )
			- ( f.expect_row_labels ? 1 : 0 );

		if( (ec = feature_alloc_encode_state( &f ) ) == 0 ) {

			if( (ec = feature_encode( line, &f, &d ) ) == 0 ) {

				_emit( line, &f, &d, stdout );
				exit_status = EXIT_SUCCESS;

			} else
				warnx( "error %d: feature_encode", ec );

			feature_free_encode_state( &f );

		} else
			warnx( "error %d: feature_alloc_encode_state", ec );

		free( line );

	} else
		warnx( "no input on stdin" );

	return exit_status;
}

#endif

