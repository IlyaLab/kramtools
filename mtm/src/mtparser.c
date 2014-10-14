
/**
  * See README.rst for format details.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <err.h>
#include <assert.h>

#ifdef HAVE_MD5
#include "md5/md5.h"
#ifndef MD5_DIGEST_LENGTH
#define MD5_DIGEST_LENGTH (16)
#endif
#endif

#include "fnv/fnv.h"

#include "syspage.h"
#include "strset.h"
#include "toktype.h"
#include "mtmatrix.h"
#include "mtsclass.h"
#include "mterror.h"

extern int mtm_sclass_by_prefix( const char *token );
extern int cardinality(
		const unsigned int *buf, int len, int largest_of_interest, const unsigned int NA );

static const char *_BUG
	= "shouldn't be here: %s:%d";

const char *MTMLIB = "libmtm";

/***************************************************************************
  * Input parameters
  */

const char *mtm_default_NA_regex = "^[Nn][Aa][Nn]?$";

#define C_STR_TERMINATOR ('\0')
#define LINE_DELIMITER   ('\n')

/**
  * These are redefinable by defining environment vars:
  * MTM_COMMENT_FLAG
  * MTM_SEPARATOR_CHAR
  */
static const char *ENVVAR_SEPARATOR = "MTM_SEPARATOR_CHAR";
static char SEPARATOR      = '\t';
static const char *ENVVAR_COMMENT   = "MTM_COMMENT_CHAR";
static char COMMENT_PREFIX = '#';

/***************************************************************************
  * Parse result/output parameters
  */


/**
  * This struct contains -all- the temporary buffers, caches, etc. required
  * for parsing...localized in one struct for (de)allocation as a unit.
  */
enum Cache {
	MATRIX_CACHE,     // Order matters! ...because _alloc_scratch only
	DESCRIPTOR_CACHE, // allocates the first 2 caches when row names
	ROWID_CACHE,      // are not going to be used.
	ROWMAP_CACHE,
	CACHE_COUNT
};

struct scratch {

	int data_column_count;
	int max_categories;

	/**
	  * The same buffer is used for both encoding types (float/int) since
	  * they're the same size and never used simultaneously.
	  * This union merely obviates ugly casting.
	  */
	union {
		mtm_fp_t  *num;
		mtm_int_t *cat;
	} buf;

	/**
	  * The string set used to map categorical features to integer values.
	  */
	void *set;

	/**
	  * tmp files to contain variable (and unpredictably!)-sized arrays.
	  */

	FILE *cache[ CACHE_COUNT ];
};

static const char *BLOCK_NAMES[] = {
	"    matrix",
	"descriptor",
	" row names",
	"   row map"
};

static void _free_scratch( struct scratch *s ) {

	// Must be tolerant of -incomplete- initialization!

	if( s->set )
		szs_destroy( s->set );

	for(int i = 0; i < CACHE_COUNT; i++ ) {
		if( s->cache[i] )
			fclose( s->cache[i] );
	}

	if( s->buf.num )
		free( s->buf.num );

	memset( s, 0, sizeof(struct scratch) );
}


/**
  * This is atomic: if it doesn not entirely succeed, it deallocates/
  * frees whatever -was- allocated and returns CLEANLY.
  */
static int _alloc_scratch( struct scratch *s,
		int max_allowed_categories,
		bool allocate_row_id_caches,
		int verbosity ) {

	const int cache_count
		= allocate_row_id_caches
		? CACHE_COUNT
		: 2;

	int i;

	assert( sizeof(mtm_int_t) == sizeof(mtm_fp_t) );
	assert( s->data_column_count > 1 );

	s->max_categories = max_allowed_categories;

	// Buffers

	s->buf.num = calloc( s->data_column_count, sizeof(mtm_fp_t) );
	if( NULL == s->buf.num ) {
		if( verbosity > 0 )
			warn( "%s: allocating buffer: %s:%d\n", MTMLIB, __FILE__, __LINE__ );
		goto failure;
	}

	// Temp file caches

	for(i = 0; i < cache_count; i++ ) {
		s->cache[i] = tmpfile();
		if( s->cache[i] == NULL )
			break;
	}
	if( i < cache_count ) {
		if( verbosity > 0 )
			warn( "%s: creating tmp file cache: %s:%d\n", MTMLIB, __FILE__, __LINE__ );
		goto failure;
	}

	// Hash table for counting categorical variables' categories.
	// Failure can only be due to out-of-memory condition.

	s->set = szs_create( max_allowed_categories,
			0, // strdup not necessary
			fnv_32_str, FNV1_32_INIT );
	if( s->set == NULL ) {
		if( verbosity > 0 )
			warn( "%s: creating string hash table: %s:%d\n", MTMLIB, __FILE__, __LINE__ );
		goto failure;
	}

	return 0;

failure:

	_free_scratch( s );
	return MTM_E_IO;
}


/**
  * Convert a string containing some number of fields into binary form,
  * possibly inferring the type of data in the process.
  *
  * Strings are mapped in the order they first appear in the input (so the
  * first string is ALWAYS mapped to 0, etc).
  */
static int _parseLine( char *line,
		struct scratch *s,
		int stat_class,
	   	int verbosity,
		struct mtm_descriptor *d ) {

	int field_type
		= field_type_from_stat_class( stat_class );
	const bool FIELD_TYPE_INFERRED
		= (__builtin_popcount( field_type ) != 1 );
		// ...either MTM_FIELD_TYPE_UNK or multiple bits.

	const char *token;
	char *pc = line;

	bool eol = false;
	unsigned int count = 0;
	int missing_value_count = 0;

	/**
	  * This union is used to quickly detect whether or not
	  */
	union {
		mtm_fp_t  f;
		mtm_int_t i;
	} last_value_read;

	assert( szs_count( s->set ) == 0 );

	while( count < s->data_column_count && ! eol ) {

		char *endpt = ""; // ESSENTIAL initialization.
		token = pc;

		while( *pc && *pc != SEPARATOR )
			pc++;

		if( token == pc /* we didn't move */ ) {

			++missing_value_count;
			s->buf.cat[count++] = NAN_AS_UINT;

			if( *pc++ )
				continue;
			else
				break; // ...because we're already ON eol.
		}

		// We have a non-empty token.
		// Before NUL-terminating see whether it's also the end-of-line.

		if( *pc == C_STR_TERMINATOR )
			eol = true; // We HAVE last field now; preclude loop reentry!
		else
			*pc++ = C_STR_TERMINATOR;

		// The "missing data" marker is common to all rows of the matrix
		// so its detection can precede (and preclude) type-dependent ops.

		if( toktype_is_na_marker(token) ) {
			++missing_value_count;
			s->buf.cat[count++] = NAN_AS_UINT;
			continue;
		}

		// For every line we will -not- know the field_type at this point
		// for at most one non-missing token. In other words, almost every
		// time we reach this point field_type is known. Use of the goto thus
		// obviates a constitutive conditional...justifying an ugly goto.

retry:
		switch( field_type ) {

		case MTM_FIELD_TYPE_FLT:
			s->buf.num[ count++ ]
				= last_value_read.f
				= strtof( token, &endpt );
			break;

		case MTM_FIELD_TYPE_INT:
			s->buf.cat[ count++ ]
				= last_value_read.i
				= strtol( token, &endpt, 0 );
			if( ! ( last_value_read.i < NAN_AS_UINT ) ) {
				return MTM_E_LIMITS;
			}
			break;

		case MTM_FIELD_TYPE_STR:
			if( szs_insert( s->set, token, &(last_value_read.i) ) == SZS_TABLE_FULL )
				return MTM_E_LIMITS;
			// The hashset will only return SZS_TABLE_FULL when genuinely
			// full, but, if s->max_categories was not itself a power of 2,
			// s->set can hold > s->max_categories. Thus, another check...
			if( szs_count( s->set ) > s->max_categories )
				return MTM_E_LIMITS;
			s->buf.cat[ count++ ] = last_value_read.i;
			break;

		default: // because either _UNKNOWN or multiple bits were set.
			// This will insure that boolean data represented by {0,1} is
			// parsed as integral data, not strings, and thus has its
			// implicit ordering preserved
			field_type = toktype_infer( token, NULL ); // can't fail.
#ifdef _DEBUG
			if( verbosity > 0 ) {
				fprintf( stderr, "%s: inferred <%s> syntax type from \"%s\"\n",
					MTMLIB, toktype_name[field_type], token );
			}
#endif
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

			if( FIELD_TYPE_INFERRED  ) {

				if( (field_type == MTM_FIELD_TYPE_INT)
					&& (MTM_FIELD_TYPE_FLT == toktype_infer( token, NULL ) ) ) {
#ifdef _DEBUG
					fputs( "promoting integral line to float\n", stderr );
#endif
					field_type = MTM_FIELD_TYPE_FLT;
					--count;
					// Convert all earlier values to float...
					for(int i = 0; i < count; i++ ) {
						if( NAN_AS_UINT != s->buf.cat[ i ] ) {
							s->buf.num[i] = (float)s->buf.cat[i]; // in place!
						}
					}
					// ...and reparse current token.
					s->buf.num[ count++ ]
						= strtof( token, &endpt );
				} else
				if( (field_type != MTM_FIELD_TYPE_STR)
					&& (MTM_FIELD_TYPE_STR == toktype_infer( token, NULL ) ) ) {
#ifdef _DEBUG
					fputs( "revising non-string line to string\n", stderr );
#endif
					assert( szs_count( s->set ) == 0 );
					field_type = MTM_FIELD_TYPE_STR;

					// Reparse line up to and including current token.
					// Above inference of eol still holds, and, importantly,
					// we can count on NUL-termination of ALL relevant
					// tokens.

					pc = line;
					count = 0;
					do {
						if( szs_insert( s->set, pc, &(last_value_read.i) ) == SZS_TABLE_FULL )
							return MTM_E_LIMITS;
						// See comment re: szs_insert elsewhere.
						if( szs_count( s->set ) > s->max_categories )
							return MTM_E_LIMITS;
						s->buf.cat[ count++ ] = last_value_read.i;
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

	if( count < s->data_column_count ) {
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
	// The first two are trivial to check for...

	if( count - missing_value_count < 2 )
		d->constant = 1;

	// ...so we only check for the 3rd if necessary (below).

	switch( field_type ) {

	case MTM_FIELD_TYPE_FLT:
		if( ! d->constant &&
				cardinality( s->buf.cat, count, 2, NAN_AS_UINT ) < 2 ) {
			d->constant = 1;
		}
		break;

	case MTM_FIELD_TYPE_STR:
		d->integral = 1;
		d->categories = szs_count( s->set );
		if( d->categories < 2 )
			d->constant = 1;
		szs_clear( s->set );
		break;

	case MTM_FIELD_TYPE_INT:
		d->integral = 1;
		// Categorical (esp. BOOLEAN categorical) features may have been encoded
		// directly in integral values instead of strings. Detect that here...
		if( ! d->constant ) {
			const int card
				= cardinality( s->buf.cat, count, s->max_categories, NAN_AS_UINT );
			// If the cardinality is small enough (as defined by input parameters)
			// it's treatable as a categorical variable whether it -is- or not...
			// TODO: Ordinal is not fully supported yet; rows containing integers
			// that the caller did not explicitly specify were continuous are
			// assumed to be categorical, so (currently) the cardinality of the
			// set of values must also respect the category limit. Eventually,
			// this may not be true.
			if( card <= s->max_categories )
				d->categories = card;
			else
				return MTM_E_LIMITS;
		}
		// If the feature is encoded as integers and more than the maximum
		// allowed categories exist, it is implicitl ordinal(?)
		break;

	case MTM_FIELD_TYPE_UNK:
	default:
		if( missing_value_count < count )
			errx( -1, _BUG, __FILE__, __LINE__ );
		// ...but if whole line was empty we may, indeed, not know the type!
	}

	return 0;
}


/**
  * Using separated fields, column counting is trivial: there will
  * AlWAYS be one more column than there are separators.
  *
  * If the line contains only separators, all fields are empty!
  *
  * Ihis, in particular, allows for the (R-style) of an empty
  * left-justified field in the header line.
  */
static int _count_columns( const char *pc ) {

	int n = 0;
	while( *pc ) if( *pc++ == SEPARATOR ) n++;
	return ++n; // ..."++" for the implicit last token.
}


/**
  * The parser allocates the matrix as one large memory blob, so
  * this destroyer is particularly simple.
  */
static void _free_matrix( struct mtm_matrix *m ) {
	if( m ) {
		if( m->storage )
			free( m->storage );
		m->storage = NULL;
	}
}


/**
  * Parse a text matrix satisfying the format description (elsewhere)
  * leaving the results in RAM and filling out the struct mtm_matrix with
  * the locations of the various parts of the parse results.
  *
  * Parsing means:
  * 1. the file is converted to binary form
  * 2. the format is implicitly validated (error-free return)
  * 3. univariate degeneracies in data rows are detected and characterized.
  *
  * This function can leave its results in one (or both) of two forms:
  * 1. memory resident, contained in a struct mtm_matrix.
  * 2. file resident in a binary file.
  *
  * The <m> is non-NULL a memory-resident result is created.
  * If <outout> is non-NULL the result is written to file.
  */
int mtm_parse( FILE *input,
		unsigned int flags,
		const char *missing_data_regex,
		int max_allowed_categories,
		int (*infer_stat_class)(const char *),
		struct mtm_matrix *m ) {

	const int verbosity
		= MTM_VERBOSITY_MASK & flags;
	const bool EXPECT_ROW_NAMES
		= ( flags & MTM_MATRIX_HAS_ROW_NAMES) != 0;
	const bool EXPECT_HEADER
		= ( flags & MTM_MATRIX_HAS_HEADER )   != 0;
	const bool PRESERVE_ROWNAMES
		= EXPECT_ROW_NAMES && ( (flags & MTM_DISCARD_ROW_NAMES) == 0 );

	int    lnum = 0;
	int    fnum = 0;

	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;

	int econd = 0;
	struct scratch s;
	struct mtm_descriptor d;

#ifdef HAVE_MD5
	md5_state_t hashstate;
	md5_byte_t checksum[ MD5_DIGEST_LENGTH ];
	md5_init( &hash_state );
#endif
	memset( &s, 0, sizeof(struct scratch));

	/**
	  * Possibly redefine SEPARATOR and COMMENT flags.
	  */

	if( getenv(ENVVAR_COMMENT) )
		COMMENT_PREFIX = getenv(ENVVAR_COMMENT)[0];
	if( getenv(ENVVAR_SEPARATOR) )
		SEPARATOR = getenv(ENVVAR_SEPARATOR)[0];

	econd = toktype_init(
		missing_data_regex ? missing_data_regex : mtm_default_NA_regex );
	if( econd ) {
		return MTM_E_INIT_REGEX;
	}

#ifdef _DEBUG
	if( verbosity > 0 ) {
		int n = 0;
		if( EXPECT_ROW_NAMES ) n += fprintf( stderr, "%s: expect row names ", MTMLIB );
		if( EXPECT_HEADER    ) n += fprintf( stderr, "%s: expect column names ", MTMLIB );
		if( n > 0 ) fputc( '\n', stderr );
	}
#endif

	/**
	  * Basic parse strategy is line-oriented:
	  * 1) read one line at a time
	  * 2) parse and process it (AFTER checksumming the unadulterated line)
	  * 3) store results in temporary file "caches"
	  * 4) after file is consumed, reconstitute the caches into RAM
	  */
	while( ( llen = getline( &line, &blen, input ) ) > 0 ) {

		int stat_class = MTM_STATCLASS_UNKNOWN;

		char *pc = line;

		++lnum; // ...at beginning of loop for 1-based line reporting.

		assert( line[llen] == C_STR_TERMINATOR ); // getline well-behaved.

		// Checksum before ANY changes happen to the read bytes...

#ifdef HAVE_MD5
		if( hashstate ) md5_append( hashstate, (md5_byte_t*)line, llen );
#endif
		// Now replace newline, if present, with NUL terminator to
		// GUARANTEE NUL termination for all subsequent code.
		// (Last line of file might not have a newline.)

		if( line[llen-1] == LINE_DELIMITER ) line[--llen] = C_STR_TERMINATOR;

		// Skip empty lines and comments.

		if( llen == 0 || line[0] == COMMENT_PREFIX )
			continue;

		// We have a "real" line. If it's the first infer column count.

		if( s.data_column_count == 0 ) {

			s.data_column_count
				= _count_columns( line ) - ( EXPECT_ROW_NAMES ? 1 : 0 );

			if( s.data_column_count < 2 /* absolute minimum sensible */ ) {
				econd = MTM_E_FORMAT_MATRIX;
				break;
			}
			if( _alloc_scratch( &s, max_allowed_categories, PRESERVE_ROWNAMES, verbosity ) ) {
				econd = MTM_E_SYS;
				break;
			}
			if( EXPECT_HEADER /* i.e. 1st "real" line is header */ )
				continue; // no further processing on this line.
		}

		if( EXPECT_ROW_NAMES ) {

			// Scan past the row identifier.

			while( *pc != SEPARATOR && /* to be safe */ *pc != LINE_DELIMITER ) pc++;

			// ...and make sure there was more than just a row id on the line!

			if( SEPARATOR != *pc ) {
				econd = MTM_E_FORMAT_MATRIX; // non-empty line with exactly one field; that's wrong!
				break;
			}

			*pc++ = C_STR_TERMINATOR;

			if( infer_stat_class )
				stat_class = infer_stat_class( line );

			if( PRESERVE_ROWNAMES ) {
				const struct mtm_row_id srn = {
					fnum,
					(const char*)ftell(s.cache[ROWID_CACHE])
				};
				const int NCHAR
					= pc-line; // include the NUL terminator!
				if( fwrite( line, sizeof(char), NCHAR, s.cache[ROWID_CACHE] )
						!= NCHAR
					||
					fwrite( &srn, sizeof(struct mtm_row_id), 1, s.cache[ROWMAP_CACHE] )
						!= 1 ) {
					econd = MTM_E_IO;
					break;
				}
			}
		}

		memset( &d, 0, sizeof(d) );

		if( ( econd = _parseLine( pc, &s, stat_class, verbosity, &d ) ) ) {
			warnx( "%s: aborting parsing at input line %d", __FILE__, lnum );
			break;
		}

		/**
		  * Finally write out the converted row and its descriptor.
		  */

		if( fwrite( s.buf.cat, sizeof(mtm_int_t), s.data_column_count, s.cache[MATRIX_CACHE] )
				!= s.data_column_count ) {
			econd = MTM_E_IO;
			break;
		}
		if( fwrite( &d, sizeof(struct mtm_descriptor), 1, s.cache[DESCRIPTOR_CACHE] )
				!= 1 ) {
			econd = MTM_E_IO;
			break;
		}

		fnum += 1;
	}

	if( line )
		free( line );

#ifdef HAVE_MD5
	md5_finish( &hash_state, checksum );
#endif

	/**
	  * Reconstitute the matrix in RAM from the tmpfile caches.
	  */

	if( m ) {

		size_t    sizeof_part[4];
		size_t pa_sizeof_part[4];
		char *blob;

		memset(    sizeof_part, 0, sizeof(   sizeof_part) );
		memset( pa_sizeof_part, 0, sizeof(pa_sizeof_part) );
		memset(              m, 0, sizeof(struct mtm_matrix) );

		// Determine the file sizes of each cache and total size

		for(int i = 0; i < CACHE_COUNT; i++ ) {
			if( s.cache[i] ) {
				sizeof_part[i]    = ftell( s.cache[i] );
				pa_sizeof_part[i] = page_aligned_ceiling( sizeof_part[i] );
				m->size          += pa_sizeof_part[i];
#ifdef _DEBUG
				if( verbosity > 0 ) {
					fprintf( stderr, "%s: cache %d: %ld bytes on disk, %ld bytes in RAM\n",
						MTMLIB, i, sizeof_part[i], pa_sizeof_part[i] );
				}
#endif
			}
		}

#ifdef _DEBUG
		if( verbosity > 0 ) {
			fprintf( stderr, "%s: %d row X %d column matrix, sizeof memory-resident rep is %ld\n",
				MTMLIB, fnum, s.data_column_count, m->size );
		}
#endif

		if( posix_memalign( (void**)&blob, RT_PAGE_SIZE, m->size ) == 0 ) {

			int i;
			char *ptr = blob;
			for(i = 0; i < CACHE_COUNT; i++ ) {
				if( s.cache[i] ) { // ignore ununsed caches

					// Load the cache back into RAM...

					if( fflush( s.cache[i] ) ||
						fseek( s.cache[i], 0, SEEK_SET ) ||
						fread( ptr, sizeof_part[i], 1, s.cache[i] ) != 1 ) {
						warnx( "%s: failed reading cache %d", __FILE__, i );
						break;
					}
#ifdef _DEBUG
					if( verbosity > 2 )
						fprintf( stderr, "%s: %s <= %p (%ld bytes)\n",
							MTMLIB, BLOCK_NAMES[i], ptr, sizeof_part[i] );
#endif
					// ...and update the matrix pointer.

					switch( i ) {
					case MATRIX_CACHE:
						m->data      = (mtm_int_t *)ptr;
						break;
					case DESCRIPTOR_CACHE:
						m->prop      = (struct mtm_descriptor *)ptr;
						break;
					case ROWID_CACHE:
						m->row_names = (const char *)ptr;
						break;
					case ROWMAP_CACHE:
						m->row_map   = (struct mtm_row_id *)ptr;
						break;
					default:
						errx( -1, _BUG, __FILE__, __LINE__ );
					}
					ptr += pa_sizeof_part[i];
				}
			}
			if( i < CACHE_COUNT ) {
				free( blob );
				blob = NULL;
			} else {
				m->destroy = _free_matrix;
				m->storage = blob;
			}
		}

		m->rows    = fnum;
		m->columns = s.data_column_count;

		if( m->row_names != NULL && m->row_map != NULL )
			mtm_resolve_rownames( m, (signed long)m->row_names );

		// Rows' order is currently the same as in the input. This might be
		// lexigraphically sorted order, but we don't *know* that, so...

		m->lexigraphic_order = false;
#ifdef HAVE_MD5
		for(int i = 0; i < MD5_DIGEST_LENGTH; i++ )
			sprintf( m->md5 + 2*i, "%02x", checksum[i] );
		m->md5[ MD5_DIGEST_LENGTH*2 ] = 0;
#endif
	}

	_free_scratch( &s );

	return econd;
}


#ifdef _UNIT_TEST_PARSE
#error "mtproc implemented by main.c is the unit test of this module"
#endif

