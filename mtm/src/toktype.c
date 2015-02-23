
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <regex.h>
#include <alloca.h>
#include <assert.h>

#include "toktype.h"
#include "mtsclass.h" // for field_type_x definitions.

extern const char *MTMLIB;

/**
  * Notice that these are strict patterns: 
  * 1. octal is required to begin with a '0', 
  * 2. decimal is required to begin with a digit other than 0, and
  * 3. hexadecimal is required to begin with 0x.
  *
  * Unfortunately sometimes decimal is front padded with zeros, too.
  */
static const char *OCT_PATTERN 
	= "^0[0-7]+$";

static const char *DEC_PATTERN 
	= "^(0|[1-9][0-9]*)$";

static const char *HEX_PATTERN 
	= "^0x[0-9a-f]+$"; // case insensitivity added by regcomp

static const char *FP_PATTERN 
	= "^([-+]?(0\\.|([1-9][0-9]*)?\\.?)[0-9]*(e[-+]?[0-9]+)?|nan|inf)$";
//	= "^([-+]?[0-9]*\\.?[0-9]+(e[-+]?[0-9]+)?|nan|inf)$";

static regex_t _rx_na;
static regex_t _rx_oct;
static regex_t _rx_dec;
static regex_t _rx_hex;
static regex_t _rx_fp;

static bool _is_valid_octal( const char *sz ) {
	return regexec( &_rx_oct, sz, 0, NULL, 0 ) == 0;
}

static bool _is_valid_decimal( const char *sz ) {
	return regexec( &_rx_dec, sz, 0, NULL, 0 ) == 0;
}

static bool _is_valid_hexadecimal( const char *sz ) {
	return regexec( &_rx_hex, sz, 0, NULL, 0 ) == 0;
}

static bool _is_valid_fp( const char *sz ) {
	return regexec( &_rx_fp,  sz, 0, NULL, 0 ) == 0;
}


static void _free() {
	regfree( &_rx_fp  );
	regfree( &_rx_hex );
	regfree( &_rx_dec );
	regfree( &_rx_oct );
	regfree( &_rx_na  );
}

const char *toktype_name[] = {
	"unknown",
	"string",
	"integral",
	"integral,string",
	"floating",
	"floating,string",
	"floating,integral",
	"floating,string,integral",
};

int toktype_init( const char *na_expression ) {

	/**
	  * Only the first might fail since it's user-supplied.
	  * The remainder are hardcoded and verified during development.
	  */

	int regerr
		= regcomp( &_rx_na, na_expression, REG_EXTENDED | REG_NOSUB );
	if( regerr ) {
		const size_t required 
			= regerror( regerr, &_rx_na, NULL, 0 );
		char *buf 
			= alloca( required );
		regerror( regerr, &_rx_na, buf, required );
		fprintf( stderr, "%s: regex error: %s\n", MTMLIB, buf );
		return -1;
	}

	regerr = regcomp( &_rx_oct, OCT_PATTERN, REG_EXTENDED | REG_NOSUB );
	assert( regerr == 0 );

	regerr = regcomp( &_rx_dec, DEC_PATTERN, REG_EXTENDED | REG_NOSUB );
	assert( regerr == 0 );

	regerr = regcomp( &_rx_hex, HEX_PATTERN, REG_EXTENDED | REG_NOSUB | REG_ICASE );
	assert( regerr == 0 );

	regerr = regcomp( &_rx_fp,   FP_PATTERN, REG_EXTENDED | REG_NOSUB | REG_ICASE );
	assert( regerr == 0 );

	atexit( _free );

	return 0;
}


bool toktype_is_na_marker( const char *sz ) {
	return regexec( &_rx_na, sz, 0, NULL, 0 ) == 0;
}


/**
  * This function must not fail. Specifically, it must return EXACTLY ONE
  * of the non-zero MTM_FIELD_TYPE_x constants.
  * "Narrowness" refers to interpretation. Specifically, everything that
  * is interpretable as an integer is interpretable as floating-point (even
  * though there are large integers that cannot be represented as floats!).
  * Obviously, everything is interpretable as a string, so that's the catch-
  * all case.
  */
int toktype_infer_narrowest_type( const char *sz, unsigned int *base ) {

	// This function should not be called on missing data markers. 
	// Such markers are, by design, common for all rows of a matrix; 
	// they cannot imply anything about the data type.

	assert( ! toktype_is_na_marker( sz ) );

	if( _is_valid_octal( sz ) ) {
		if( base ) *base =  8;
		return MTM_FIELD_TYPE_INT;
	} else
	if( _is_valid_decimal( sz ) ) {
		if( base ) *base = 10;
		return MTM_FIELD_TYPE_INT;
	} else
	if( _is_valid_hexadecimal( sz ) ) {
		if( base ) *base = 16;
		return MTM_FIELD_TYPE_INT;
	} else
	if( _is_valid_fp( sz ) ) {
		return MTM_FIELD_TYPE_FLT;
	}

	return MTM_FIELD_TYPE_STR;
}

#ifdef _UNIT_TEST_TOKTYPE

#include <err.h>

int main( int argc, char *argv[] ) {

	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;
	if( argc > 1 && toktype_init( argv[1] ) == 0 ) {
		int tt;
		while( ( llen = getline( &line, &blen, stdin ) ) > 0 ) {
			unsigned int base;
			if( line[llen-1] == '\n' ) line[--llen] = '\0';
			tt = toktype_infer( line, &base );
			if( tt == MTM_FIELD_TYPE_INT )
				printf( "%s (%d)\n", toktype_name[ tt ], base );
			else
				puts( toktype_name[ tt ] );
		}
		if( line )
			free( line );
	} else
		errx( -1, "failed initializing" );
	return 0;
}
#endif

