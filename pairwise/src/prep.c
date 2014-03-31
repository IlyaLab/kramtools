
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "strset.h"

/**
  * This file is (almost) a drop-in compiled replacement for prep.py.
  * "Almost" because if a file is being pre-processed just-in-time then
  * some of the persistent structure created by prep.py need not be
  * created; this is not creating a preprocessed file that will be
  * used repeatedly.
  *
  * Definition of valid input format:
  * 1. File contains '\t'-delimited fields and '\n'-delimited lines.
  * 2. No other character satisfying isspace is allowed anywhere in
  *    the file. In particular: 
  *    a. "normal" (' ') spaces are illegal. With minor changes they 
  *       -could- be allowed but they're such a Bad Idea in general
  *       that it's just not worthwhile.
  *    b. '\r' is also disallowed. In other words: Unix line endings.
  * 3. Tabs are expected as delimiters, not merely separators. What's the
  *    difference? Adjacent tabs are NOT coalesced into a single separating 
  *    gap. For example, three adjacent tabs imply two empty fields.
  *
  * It converts to binary and preprocesses \r?\n-delimited lines of 
  * \t-delimited fields.
  */
#define MAX_ALLOWED_CATEGORY_COUNT (32)
#define C_STR_TERMINATOR ('\0')
#define LINE_DELIMITER   ('\n')
#define FIELD_DELIMITER  ('\t')
#define NAN_AS_UINT (0x7FC00000)

#define SIZEOF_ELEMENT  (sizeof(float))
/**
  * ONE-PASS parsing, validation, and conversion to binary is required.
  */


/**
  * <N> Number of columns to expect.
  * As noted elsewhere, this assumes TABs are delimiters, so
  * 1. two consecutive TABs constitute an empty field and 
  * 2. will result in failure to read N fields.
  *
  * WARNING: NO PART OF THIS FUNCTION MAY RELY ON NUL-TERMINATION OF LINE.
  */
static int _processNumeric( char *line, uint32_t *tag, 
		const int EXPECTED_FIELD_COUNT, float *buf ) {

	const char *pc = line;
	int count  = 0;
	int missing_value_count = 0;
	float last_normal_value;

	while( count < EXPECTED_FIELD_COUNT ) {

		assert( isgraph(*pc) ); // pc should be on -start- of token

		// Following mutually exclusive possibilities effectively 
		// disallow "nan" from appearing in matrix.

		if( isalpha(*pc) ) {
			if( strncasecmp( pc, "NA", 2 ) == 0 ) {
				++missing_value_count;
				buf[count++] = nanf("nan");
				pc += 2;
			} else
				return -1;
		} else {
			last_normal_value
				= buf[count++] 
				= strtof( pc, (char **)(&pc) );
			if( ERANGE == errno )
				return -1;
		}

		// pc should be on 1st char after token, and it MUST
		// be either a line or field delimiter.

		if( LINE_DELIMITER  == *pc )
			break;
		else
		if( FIELD_DELIMITER != *pc++ )
			return -1;
	}

	if( count < EXPECTED_FIELD_COUNT ) {
		return -2;
	}

	// Now check for the various kinds of degeneracy:
	// 1. Too few (<2) values present (non-missing)
	// 2. |{x:x is non-missing}| == 1, i.e. a "constant" vector.

	if( count - missing_value_count < 2 ) {
		*tag = 0x00010000;
	} else {
		int i;
		// Insure there are at least two distinct values.
		for(i = 0; i < EXPECTED_FIELD_COUNT; i++ ) {
			const float F = buf[i];
			if( ( ! isnan(F) ) && last_normal_value != F )
				break;
		}
		if( i == EXPECTED_FIELD_COUNT )
			*tag = 0x00010000;
	}

	*tag |= (0x0000FFFF & missing_value_count);

	return 0;
}


/**
  *
  * WARNING: NO PART OF THIS FUNCTION MAY RELY ON NUL-TERMINATION OF LINE.
  */
static int _processCategorical( char *line, uint32_t *tag, 
		const int EXPECTED_FIELD_COUNT, uint32_t *buf, void *set ) {

	bool final = false;
	int last_normal_value, count = 0;
	unsigned int h;
	const char *token;
	char *pc = line;
	int missing_value_count = 0;

	assert( szs_count( set ) == 0 );

	while( count < EXPECTED_FIELD_COUNT && ! final ) {

		assert( isgraph(*pc) ); // pc should be on -start- of token

		token = pc;

		// Find current token's end.

		while( *pc != FIELD_DELIMITER && *pc != LINE_DELIMITER )
			pc++;

		// Before NUL-terminating see whether it's also the end-of-line.

		final = ( LINE_DELIMITER == *pc );
		*pc++ = '\0';

		if( strcasecmp( "NA", token ) == 0 ) {
			++missing_value_count;
			buf[count++] = NAN_AS_UINT;
		} else {
			if( szs_insert( set, token, &last_normal_value ) == SZS_TABLE_FULL ) {
				puts( "full table" );
				return -1;
			}
			buf[ count++ ] = last_normal_value;
		}
	}

	if( count < EXPECTED_FIELD_COUNT ) {
		return -2;
	}

	*tag |= (0x0000FFFF & missing_value_count);

	return 0;
}


/**
  * This counter allows the (R-style) of an empty first (upper-left) field.
  * That is, the matrix' first (non-comment, non-empty) line may begin with
  * a tab character.
  * WARNING: NO PART OF THIS FUNCTION MAY RELY ON NUL-TERMINATION OF LINE.
  */
static int _count_columns( const char *pc ) {

	int n = 0;
	while( true ) {
		const int C = *pc++;
		if( FIELD_DELIMITER == C )
			n++;
		else
		if( LINE_DELIMITER  == C ) {
			n++;
			break;
		}
	}
	return n;
}

/**
  * This struct contains -all- the temporary heap allocations required
  * for preprocessing, within preprocess(). ...for (de)allocation as a unit.
  */
struct scratch {

	/**
	  * This does -not- include the column of row identifiers.
	  * Note that both input and the resulting matrix have one
	  * more column than the count of data columns: the input
	  * matrix' first column is row identifiers, and the output
	  * matrix' first column is flags.
	  */
	int data_column_count;

	union {
		float    *num;
		uint32_t *cat;
	} buf;

	/**
	  * The string set used to map categorical features to integer values.
	  */
	void *set;

	/**
	  * The resulting binary matrix.
	  */
	FILE *matrix;
};

static void _free_scratch( struct scratch *s ) {

	szs_destroy( s->set );

	fclose( s->matrix );

	if( s->buf.num )
		free( s->buf.num );

	memset( s, 0, sizeof(struct scratch) );
}


/**
  * This is transactional: entirely succeeds or fails WITH cleanup.
  */
static int _alloc_scratch( struct scratch *s ) {

	assert( sizeof(uint32_t) == sizeof(float) );
	assert( s->data_column_count > 1 );

	// Buffers

	s->buf.num = calloc( s->data_column_count + 1, sizeof(float) );
	if( NULL == s->buf.num ) { // very unlikely
		_free_scratch( s );
		fprintf( stderr, "out of memory: %s:%d\n",
				__FILE__, __LINE__ );
		return -1;
	}

	// Accumulator

	s->matrix = tmpfile();
	if( s->matrix == NULL ) {
		_free_scratch( s );
		fprintf( stderr, "failed creating temporary file: %s:%d\n",
				__FILE__, __LINE__ );
		return -1;
	}

	// Hash table for counting categorical variables' categories.

	s->set = szs_create( MAX_ALLOWED_CATEGORY_COUNT, 0 /* strdup not necessary */ );
	if( s->set == NULL ) {
		_free_scratch( s );
		fprintf( stderr, "failed creating hash table: %s:%d\n",
				__FILE__, __LINE__ );
		return -1;
	}

	return 0;
}


int preprocess( FILE *fp, bool expect_header ) {

	int    lnum = 0;
	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;
	int error;
	int binary_column_count = 0;

	struct scratch s;
	memset( &s, 0, sizeof(struct scratch));

	/**
	  * Every line should:
	  * 1. contain data, or
	  * 2. be a #-prefixed comment, or
	  * 3. be empty (meaning containing exactly 1 newline char).
	  * Exactly one line may be a header and it must be the first non-empty,
	  * non-comment line. If a header is expected, the first non-empty, non-
	  * comment line is discarded after column counting. If a header is not
	  * expected this line is assumed to contain data.
	  */

	while( ( llen = getline( &line, &blen, fp ) ) > 0 ) {

		char *pc = line;

		++lnum; // ...at beginning of loop for 1-based line reporting.

		assert( line[llen] == C_STR_TERMINATOR ); // getline well-behaved.
		assert( llen < blen ); // buffer has room for NUL-terminator. 

		// The LAST line of input might not be newline-terminated. 
		// All the _processXXX functions ASSUME newline termination.
		// Following wart GUARANTEES "\n" termination of every line 
		// AT THE EXPENSE of '\0' termination.
		// _processXXX MUST NOT assume '\0' termination.

		if( line[llen-1] != LINE_DELIMITER )
			line[llen++]  = LINE_DELIMITER; // overwriting the NUL-terminator.

		// Skip empty lines and comments.

		if( llen == 1 /* NL delimiter only */ || line[0] == '#' )
			continue;

		// We have a "real" line. Infer column count now if haven't already.

		if( s.data_column_count == 0 ) {
			s.data_column_count = _count_columns( line ) - 1;
			binary_column_count = s.data_column_count + 1;
			if( s.data_column_count < 2 /* absolute minimum sensible */ )
				return -2;
			if( _alloc_scratch( &s ) )
				return -1;
			if( expect_header /* i.e. 1st "real" line is header */ )
				continue; // no further processing on this line.
		}

		// Scan past the row identifier.

		while( *pc != FIELD_DELIMITER && /* to be safe */ *pc != LINE_DELIMITER ) pc++;

		// ...and make sure there was more than just a row id on the line!

		if( FIELD_DELIMITER != *pc )
			return -1; // non-empty line with 1 field.

		*pc++ = C_STR_TERMINATOR;

		if( ! isgraph(*pc) ) // _processXXX functions would catch this, but
			return -1;       // in the spirit of fail-early, check it here!

		switch( toupper(line[0]) ) {
		case 'N': // most common first!
			error = _processNumeric(
					pc, s.buf.cat, s.data_column_count, s.buf.num + 1 );
			break;
		case 'B':
		case 'C':
			error = _processCategorical(
					pc, s.buf.cat, s.data_column_count, s.buf.cat + 1, s.set );
			szs_clear( s.set );
			break;
		default:
			fprintf( stderr, "unexpected line prefix: \t%s on line %d", pc, lnum );
		}
		if( ! error ) {
			if( fwrite( s.buf.cat, SIZEOF_ELEMENT, binary_column_count, s.matrix ) != binary_column_count )
				return -1;
		} else
			printf( "line %d => %d\n", lnum, error );
	}

cleanup:
	if( line )
		free( line );

	printf( "sizeof tmpfile = %ld\n", ftell( s.matrix ) );

	_free_scratch( &s );

	return 0;
}

#ifdef _UNIT_TEST_PREP
int main( int argc, char *argv[] ) {

	const int EXPECT_HEADER 
		= argc > 1 && tolower(argv[1][0]) == 'h';
	FILE *fp
		= getenv("PREP_FILE")
		? fopen( getenv("PREP_FILE"), "r" )
		: stdin;
	int status
		= preprocess( fp, EXPECT_HEADER );
	fclose( fp );

	return status;
}
#endif

