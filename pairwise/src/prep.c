
/**
  * This file is (almost) a drop-in compiled replacement for prep.py.
  * "Almost" because if a file is being pre-processed just-in-time then
  * some of the persistent structure created by prep.py need not be
  * created; this is not creating a preprocessed file that will be
  * used repeatedly.
  *
  * Definition of valid input format:
  * 1. File contains '\t'-delimited fields and '\n'-delimited lines.
  * 2. Empty fields are illegal.
  * 3. No other character satisfying isspace is allowed anywhere in
  *    the file. In particular: 
  *    a. "normal" (' ') spaces are illegal. With minor changes they 
  *       -could- be allowed but they're such a Bad Idea in general
  *       that it's just not worthwhile.
  *    b. '\r' is also disallowed. In other words: Unix line endings.
  * 4. Tabs are expected as delimiters, not merely separators. What's the
  *    difference? Adjacent tabs are NOT coalesced into a single separating 
  *    gap. For example, two adjacent tabs implies an empty field and, per
  *    #2, is illegal.
  *
  * The preceding apply to ALL lines of input.
  * Additionally:
  * 5. Each line may be 
  *    a. empty (consisting of exactly one '\n' character)
  *    b. a comment if the first character is '#'
  *    c. a '\t'-delimited data or (also '\t'-delimited) header line.
  *
  * 6. All data and header lines must contain the same column count
  *    (which means the same number of '\t' characters).
  * 7. A header, if present, must be the first non-empty, non-comment
  *    line.
  * 8. Row identifiers are optional, but if present they will always
  *    be left-justified (NO SPACES) in the first column.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "limits.h"
#include "strset.h"
#include "featqual.h"
#include "preprohd.h"
#include "matrix.h"

#ifdef HAVE_CHECKSUM
#include "md5/md5.h"
#ifndef MD5_DIGEST_LENGTH
#define MD5_DIGEST_LENGTH (16)
#endif
#endif

#define C_STR_TERMINATOR ('\0')
#define LINE_DELIMITER   ('\n')
#define FIELD_DELIMITER  ('\t')

#define INDICATOR_NUMERIC     'N'
#define INDICATOR_BOOLEAN     'B'
#define INDICATOR_CATEGORICAL 'C'

#define NAN_AS_UINT (0x7FC00000)
#define SIZEOF_HEADER  (96)

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
static int _processNumeric( char *line, struct feature_qualifier *fq, 
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

	// Now check for the degeneracy. There is really only one kind of
	// degeneracy, constancy, which can manifest in three different ways:
	// 1. All values are missing, so feature is entirely "NA"
	// 2. Exactly one-value is NOT missing, so obv. there can only be one value.
	// 3. Two or more values are non-missing, but they are all equal.

	if( count - missing_value_count < 2 )
		fq->constant = 1;
	else {
		int i;
		// Insure there are at least two distinct values.
		for(i = 0; i < EXPECTED_FIELD_COUNT; i++ ) {
			const float F = buf[i];
			if( ( ! isnan(F) ) && last_normal_value != F )
				break;
		}
		if( i == EXPECTED_FIELD_COUNT )
			fq->constant = 1;
	}

	fq->missing = missing_value_count;

	return 0;
}


/**
  * Convert a line of strings, arbitrary except that they should contain 
  * NO whitespace, into an array of non-negative integers such that each
  * string is mapped to a unique integer.
  *
  * One string in the input is special: /[nN][aA]/. "NA" is mapped to the
  * floating point NAN value (0x7FC00000) which is cast to unsigned int.
  *
  * Strings are mapped in the order they first appear in the input (so the
  * first string is ALWAYS mapped to 0, etc).
  *
  * After mapping, check for 
  * WARNING: NO PART OF THIS FUNCTION MAY RELY ON NUL-TERMINATION OF LINE.
  */
static int _processCategorical( char *line, struct feature_qualifier *fq, 
		const int EXPECTED_FIELD_COUNT, unsigned int *buf, void *set, FILE *meta ) {

	bool final = false;
	unsigned int last_normal_value, count = 0;
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
				return __LINE__;
			}
			buf[ count++ ] = last_normal_value;
		}
	}

	if( count < EXPECTED_FIELD_COUNT ) {
		return __LINE__;
	}

	fq->categories = szs_count( set );
	fq->missing    = missing_value_count;

	if( count - missing_value_count < 2 || fq->categories == 1 )
		fq->constant = 1;
	else {
		int i;
		// Insure there are at least two distinct values.
		for(i = 0; i < EXPECTED_FIELD_COUNT; i++ ) {
			const unsigned int I = buf[i];
			if( (I != NAN_AS_UINT) && (I != last_normal_value) )
				break;
		}
		if( i == EXPECTED_FIELD_COUNT )
			fq->constant = 1;
	}

	return 0;
}


/**
  * This counter allows the (R-style) of an empty first (upper-left) field.
  * That is, the matrix' first (non-comment, non-empty) line may begin with
  * a tab character.
  *
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
		unsigned int *cat;
	} buf;

	struct feature_qualifier *qualifier;

	/**
	  * The string set used to map categorical features to integer values.
	  */
	void *set;

	/**
	  * The resulting binary matrix.
	  */
	FILE *matrix;
	FILE *rownames;
	FILE *rowmap;
};


static void _free_scratch( struct scratch *s ) {

	szs_destroy( s->set );

	if( s->rowmap )
		fclose( s->rowmap   );
	if( s->rownames )
		fclose( s->rownames );

	fclose( s->matrix   );

	if( s->buf.num )
		free( s->buf.num );

	memset( s, 0, sizeof(struct scratch) );
}


/**
  * This is transactional: entirely succeeds or fails WITH cleanup.
  */
static int _alloc_scratch( struct scratch *s ) {

	assert( sizeof(unsigned int) == sizeof(float) );
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
	if( s->matrix == NULL || fseek( s->matrix, sizeof(struct preprocessed_input_header), SEEK_SET ) ) {
		_free_scratch( s );
		fprintf( stderr, "failed creating (or seeking in) temporary file for matrix: %s:%d\n",
				__FILE__, __LINE__ );
		return -1;
	}

	s->rownames = tmpfile();
	if( s->rownames == NULL ) {
		_free_scratch( s );
		fprintf( stderr, "failed creating temporary file for rownames: %s:%d\n",
				__FILE__, __LINE__ );
		return -1;
	}

	s->rowmap   = tmpfile();
	if( s->rowmap == NULL ) {
		_free_scratch( s );
		fprintf( stderr, "failed creating temporary file for rowmap: %s:%d\n",
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

	// ...just to avoid lots of ugly casting warts elsewhere.

	s->qualifier = (struct feature_qualifier*)s->buf.cat;

	return 0;
}


int _write_binary( struct scratch *s, int fnum, md5_state_t *hash, FILE *output ) {

	int i;
	long offset 
		= sizeof(struct preprocessed_input_header)
		+ ftell( s->matrix );
	struct preprocessed_input_header hdr;

#ifdef HAVE_CHECKSUM
	md5_byte_t checksum[ MD5_DIGEST_LENGTH ];
	md5_finish( hash, checksum );
#endif

	memset( &hdr, 0, sizeof(hdr) );

	while( offset & 0x0FFF /* PAGE_SIZE alignment for mem mapping */ ) 
		offset += 1;

	strncpy( hdr.sig, "PAIRWISE", 8 );
	hdr.version      = 0x00020000;
	hdr.header_size  = sizeof(struct preprocessed_input_header);
	hdr.count_rows   = fnum;
	hdr.offset_names = offset;

	while( offset & 0x0FFF /* PAGE_SIZE alignment for mem mapping */ ) 
		offset += 1;

	hdr.offset_map   = offset;
	hdr.flags        = 0;

#ifdef HAVE_CHECKSUM
	for(i=0; i < MD5_DIGEST_LENGTH; i++ ) 
		sprintf( hdr.md5 + 2*i, "%02x", checksum[i] );
#endif

	fseek( s->matrix, 0, SEEK_SET );
	fwrite( &hdr, sizeof(hdr), 1, s->matrix );

	return 0;
}


#define MATRIX_HAS_HEADER    (0x00000001)
#define MATRIX_HAS_ROW_NAMES (0x00000002)
#define DISCARD_ROW_NAMES    (0x00000004)

/**
  * Preprocesses a text matrix satisfying the format description at top.
  * 
  * Preprocessing means:
  * 1. the file is converted to binary form
  * 2. the format is implicitly validated (error-free return)
  * 3. univariate degeneracies in data rows are detected and characterized.
  *
  * Originally input matrices were "TCGA format", though that no longer really applies.
  *
  * This function can leave its results in one (or both) of two forms:
  * 1. memory resident, contained in a struct matrix.
  * 2. file resident in a binary file.
  *
  * The <m> is non-NULL a memory-resident result is created.
  * If <outout> is non-NULL the result is written to file.
  */
int preprocess( FILE *input, unsigned int flags, struct matrix *m, FILE *output ) {

	const bool PRESERVE_ROWNAMES
		= true;

	const bool EXPECT_ROW_NAMES
		= ( flags & MATRIX_HAS_ROW_NAMES) != 0;
	const bool EXPECT_HEADER
		= ( flags & MATRIX_HAS_HEADER ) != 0;

	int    lnum = 0;
	int    fnum = 0;
	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;
	int errnum = 0;
	int binary_column_count = 0;

	struct scratch s;

#ifdef HAVE_CHECKSUM
	md5_state_t context;
#endif

#ifdef HAVE_CHECKSUM
	if( output ) md5_init( &context );
#endif

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

	while( ( llen = getline( &line, &blen, input ) ) > 0 ) {

		char *pc = line;

		++lnum; // ...at beginning of loop for 1-based line reporting.

		assert( line[llen] == C_STR_TERMINATOR ); // getline well-behaved.

		// Checksum before ANY changes happen to the read bytes...

#ifdef HAVE_CHECKSUM
		if( output ) md5_append( &context, (md5_byte_t*)line, llen );
#endif
		//                    !!! WARNING !!!

		// The LAST line of input might not be newline-terminated. 
		// ALL code following this point ASSUMES newline termination, so
		// the following wart GUARANTEES '\n' termination of every line 
		// AT THE EXPENSE of '\0' termination.
		// Thus, NO CODE MAY ASSUME '\0' TERMINATION.
		// (And I've favoring '\n' termination over '\0' termination because
		//  this way isspace can be used to catch trouble.)

		//                    !!! WARNING !!!

		if( line[llen-1] != LINE_DELIMITER )
			line[llen++]  = LINE_DELIMITER; // overwriting the NUL-terminator.

		// Skip empty lines and comments.

		if( llen == 1 /* NL delimiter only */ || line[0] == '#' )
			continue;

		// We have a "real" line. Infer column count now if haven't already.

		if( s.data_column_count == 0 ) {

			s.data_column_count
				= _count_columns( line ) 
				- (EXPECT_ROW_NAMES?1:0);

			binary_column_count = s.data_column_count + 1;

			if( s.data_column_count < 2 /* absolute minimum sensible */ ) {
				errnum = __LINE__;
				break;
			}
			if( _alloc_scratch( &s ) ) {
				errnum = __LINE__;
				break;
			}
			if( EXPECT_HEADER /* i.e. 1st "real" line is header */ )
				continue; // no further processing on this line.
		}

		if( EXPECT_ROW_NAMES ) {

			// Scan past the row identifier.

			while( *pc != FIELD_DELIMITER && /* to be safe */ *pc != LINE_DELIMITER ) pc++;

			// ...and make sure there was more than just a row id on the line!

			if( FIELD_DELIMITER != *pc ) {
				errnum = __LINE__; // non-empty line with exactly one field; that's wrong!
				break;
			}

			*pc++ = C_STR_TERMINATOR;

			if( PRESERVE_ROWNAMES ) {
				unsigned int pair[2] = {0,fnum};
				pair[0] = ftell(s.rownames);
				fwrite( line, sizeof(char), 
					pc-line, // include the NUL terminator! 
					s.rownames );
				fwrite( pair, sizeof(unsigned int), 2, s.rowmap );
			}
		}

		if( ! isgraph(*pc) ) { // _processXXX functions would catch this, but
			errnum = __LINE__; // Check it here in order to fail-early!
			break;
		}

		*s.buf.cat = 0; // clear the feature_qualifier, 1st array element.

		/**
		  * Parse and process according to feature class.
		  */

		switch( toupper(line[0]) ) {
		case INDICATOR_NUMERIC: // most common first!
			errnum = _processNumeric( pc, s.qualifier,
					s.data_column_count, s.buf.num + 1 ) ? __LINE__ : 0;
			break;

		case INDICATOR_BOOLEAN:
		case INDICATOR_CATEGORICAL:
			errnum = _processCategorical( pc, s.qualifier,
					s.data_column_count, s.buf.cat + 1, s.set, NULL ) ? __LINE__ : 0;
			szs_clear( s.set );
			break;

		default:
			fprintf( stderr, "unexpected line prefix: \t%s on line %d", pc, lnum );
		}

		/**
		  * Finally write out the converted row.
		  */

		if( errnum ) 
			break;
		else {
			if( fwrite( s.buf.cat, SIZEOF_ELEMENT, binary_column_count, s.matrix ) != binary_column_count ) {
				errnum = __LINE__;
				break;
			}
		}
		fnum += 1;
	}

	if( line )
		free( line );
#if 0
	printf( 
		"sizeof matrix   = %ld ( =?= %d + %d TOTAL cols x %d rows x 4 == %ld )\n"
		"sizeof rownames = %ld\n"
		"sizeof rowmap   = %ld ( =?= %d x 8 == %d) \n", 
		ftell( s.matrix ) , SIZEOF_HEADER, s.data_column_count+1, fnum, SIZEOF_HEADER + (sizeof(float)*(s.data_column_count+1)*fnum),
		ftell( s.rownames ), 
		ftell( s.rowmap ),
		fnum, fnum*8 );
#endif

	/**
	  * Store the preprocessed file if this was requested.
	  */

	if( output ) {
		_write_binary( &s, fnum, &context, output );
	}

	_free_scratch( &s );

	return errnum;
}


#ifdef _UNIT_TEST_PREP

#include <alloca.h>

int main( int argc, char *argv[] ) {

	int status;
	const unsigned int FLAGS 
		= atoi( argv[1] );
	const char *FNAME = argc > 2 ? argv[2] : NULL;
	FILE *fp = FNAME ? fopen( FNAME, "r" ) : stdin;
	struct matrix m;

	memset( &m, 0, sizeof(struct matrix) );

	assert( sizeof(struct preprocessed_input_header) == SIZEOF_HEADER );

#if ! defined(HAVE_CHECKSUM)
	if( FNAME ) {
		const size_t L
			= strlen(FNAME)+8 < 33 ? 33 : strlen(FNAME)+8;
		char *buf
			= alloca( L );
		FILE *proc;

		sprintf( buf, "md5sum %s", FNAME );
		proc = popen( buf, "r" );
		if( proc ) {
			fgets( buf, 32+1, proc );
			fprintf( stdout, "MD5: %s\n", buf );
			fclose( proc );
		} else
			fprintf( stderr, "failed: %s\n", buf );
	}
#endif

	status = preprocess( fp, FLAGS, &m, (FILE*)0x1 );

	assert( sizeof(struct feature_qualifier) == sizeof(unsigned int) );
	fclose( fp );

	return status;
}
#endif

