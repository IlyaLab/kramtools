
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>    // isnan
#include <ctype.h>   // toupper
#include <string.h>
#include <getopt.h>
#include <err.h>

#include "syspage.h"
#include "mtmatrix.h"
#include "mtheader.h"

extern int mtm_sclass_by_prefix( const char *token );

#ifdef HAVE_BINARY_PERSISTENT_FORMAT
static int _save( struct mtm_matrix *m, FILE *fp ) {

	const size_t BASE
		= page_aligned_ceiling( sizeof(struct mtm_matrix_header) );
	struct mtm_matrix_header hdr;
	memset( &hdr, 0, sizeof(hdr) );

	memcpy( hdr.sig, "MULTIMAT", 8 );
	hdr.endian  = 0x04030201;
	hdr.version = 0x00020000;
	hdr.header_size = sizeof(struct mtm_matrix_header);
	hdr.datum_size  = sizeof(mtm_int_t);
	hdr.rows    = m->rows;
	hdr.columns = m->columns;

/*
	hdr.off_data     = BASE;
	hdr.off_descrip  = BASE + (char*)m->prop - (char*)m->data;
	hdr.off_string   = hdr.off_descrip + m->row_names - (char*)m->prop;
	hdr.off_rowmap   = hdr.off_string + ;
*/
	//hdr.flags;        // none defined yet

	memcpy( hdr.md5, checksum, 32 );

	if( fwrite( &hdr, sizeof(struct mtm_matrix_header), 1, fp ) != 1 )
		return -1;
	if( fwrite( m->data, m->size, 1, fp ) != 1 )
		return -1;
	return 0;
}
#endif


static void _dump( struct mtm_matrix *m, 
		const char *label_format,
		const char *float_format, 
		FILE *fp ) {

	const char *MISSING_DATA = "NA";

	MTM_INT_T *pdata = m->data;

	// Start echoing.
	// We never echo the header; it's not preserved in parsing.

	for(int r = 0; r < m->rows; r++ ) {

		const int MISSING
			= m->prop[r].missing;
		fprintf( fp, "%s%s%c:%c:%d:%d",
			m->row_map ? m->row_map[r].string : "",
			m->row_map ? "\t"                 : "",
			MISSING < m->columns ? "FI"[ m->prop[r].integral ] : '?',
			m->prop[r].constant ? '!' : '-',
			m->prop[r].categories,
			MISSING );

		if( m->prop[r].integral ) {

			for(int c = 0; c < m->columns; c++ ) {
				MTM_INT_T I = *pdata++;
				if( I == NAN_AS_UINT )
					fprintf( fp, "\t%s", MISSING_DATA );
				else
					fprintf( fp, "\t%d", I );
			}

		} else {

			MTM_FP_T *f
				= (MTM_FP_T*)pdata;
			for(int c = 0; c < m->columns; c++ ) {
				MTM_FP_T F = *f++;
				if( isnan(F) )
					fprintf( fp, "\t%s", MISSING_DATA );
				else {
					fputc( '\t', fp );
					fprintf( fp, float_format, F );
				}
			}
			pdata = (mtm_int_t*)f;

		}
		fputc( '\n', fp );
	}
}


/**
  * Input options
  */

static bool opt_expect_rownames = true;
static bool opt_expect_header   = true;
static const char *opt_missing_marker = NULL;
static int  opt_max_categories  = 32;
static int (*opt_interpret_type)( const char *token ) = mtm_sclass_by_prefix;

/**
  * Output options
  */

#ifdef HAVE_MD5
static bool opt_include_md5     = false;
#endif
static bool opt_echo            = false;
// This is the minimum precision that random.rb generates which
// is the MOST we can assume is available...
static const char *opt_float_format = "%.1e";
static const char *opt_label_format = "%s\t%d:%d:%d:%d";
static char opt_binary_file[ FILENAME_MAX+1 ];
static int opt_verbosity = 0;

static const char *USAGE = 
"%s [ options ] < base_filenames\n"
"Options:\n"
"   --rownames | -r   do NOT expect input to have row names \n"
"   --header   | -h   do NOT expect input to havea header \n"
"   --missing  | -m   set the \"missing data\" regex [ \"%s\" ]\n"
"   --maxcats  | -k   set the maximum number of categories allowed \n"
"                     categorical variables [%d]\n"
"   --infer    | -A   infer statistical class from syntax\n"
"                     Use ISB prefix convention on row labels by default.\n"
"Output options:\n"
"  --echo      | -e   echo the parsed and validated matrix to stdout [%s]\n"
#ifdef HAVE_MD5
"  --checksum  | -c   include the (MD5) checksum of the input[%s]\n"
#endif
"  --float     | -f   a C printf format string for floating point\n"
"                     display [\"%s\"]\n"
"  --label     | -l   a C printf format string for display of the\n"
"                     row name and descriptor (see Notes). [\"%s\"]\n"
"  --output    | -o   name of file to receive binary representation\n"
"                     of the input matrix. This is written only if\n"
"                     an output name is provided.\n"
"  --verbosity | -v verbosity [%d]\n"
"Notes:\n"
#ifdef _DEBUG
"This is a _DEBUG build.\n"
#endif
"Written by rkramer@systemsbiology.org.\n";


static void _print_usage( const char *exename, FILE *fp ) {
	static const char *_T = "yes";
	static const char *_F = "no";
	fprintf( fp, USAGE, 
		exename,
		opt_missing_marker,
		opt_max_categories,

#ifdef HAVE_MD5
		opt_include_md5 ? _T : _F,
#endif
		opt_echo ? _T : _F,
		opt_float_format,
		opt_label_format,
		opt_verbosity );
}


int main( int argc, char *argv[] ) {

	unsigned int flags = 0;
	FILE *input = stdin;
	int errnum;
	const char *FNAME = NULL;
	struct mtm_matrix m;

	opt_missing_marker = mtm_default_NA_regex;
	memset( opt_binary_file, 0, sizeof(opt_binary_file) );
	memset( &m, 0, sizeof(struct mtm_matrix) );

	do {
		static const char *CHAR_OPTIONS 
#ifdef HAVE_MD5
			= "rhm:k:Aecf:l:o:v:?";
#else
			= "rhm:k:Aef:l:o:v:?";
#endif
		static struct option LONG_OPTIONS[] = {

			{"rownames",   0,0,'r'},
			{"header",     0,0,'h'},
			{"missing",    1,0,'m'},
			{"maxcats",    1,0,'k'},
			{"infer",      0,0,'A'},

			{"echo",       0,0,'e'},
#ifdef HAVE_MD5
			{"checksum",   0,0,'c'},
#endif
			{"float",      1,0,'f'},
			{"label",      1,0,'l'},
			{"output",     1,0,'o'},

			{"verbosity",  1,0,'v'},
			{ NULL,        0,0, 0 }
		};

		int opt_offset = 0;
		const int c = getopt_long( argc, argv, CHAR_OPTIONS, LONG_OPTIONS, &opt_offset );
		switch (c) {

		case 'r': opt_expect_rownames = false;        break;
		case 'h': opt_expect_header   = false;        break;
		case 'm': opt_missing_marker  = optarg;       break;
		case 'k': opt_max_categories  = atoi(optarg); break;
		case 'A': opt_interpret_type  = NULL;         break;
		case 'e': opt_echo            = true;         break;
#ifdef HAVE_MD5
		case 'c': opt_include_md5     = true;         break;
#endif
		case 'f': opt_float_format = optarg;          break;
		case 'l': opt_label_format = optarg;          break;

		case 'o':
			strcpy( opt_binary_file, optarg );
			break;

		case 'v':
			opt_verbosity = atoi( optarg );
			break;

		case '?':
			_print_usage( argv[0], stdout );
			exit(0);
			break;
		case -1: // ...signals no more options.
			break;
		default:
			printf ("error: unknown option: %c\n", c );
			exit(-1);
		}
		if( -1 == c ) break;

	} while( true );

	flags
		= ( MTM_VERBOSITY_MASK & opt_verbosity)
		| ( opt_expect_rownames ? MTM_MATRIX_HAS_ROW_NAMES : 0 )
		| ( opt_expect_header   ? MTM_MATRIX_HAS_HEADER    : 0 );

	if( optind < argc ) {
		FNAME = argv[optind++];
		input = fopen( FNAME, "r" );
		if( NULL == input ) {
			err( -1, "loading %s", FNAME );
		}
	}

	errnum = mtm_parse( input, 
		flags, 
		opt_missing_marker, 
		opt_max_categories, 
		opt_interpret_type, // may be NULL
		&m );
	fclose( input );

#ifdef HAVE_BINARY_PERSISTENT_FORMAT

	// Save the binary first (before resolving the string offsets in the 
	// row_map to string pointers).

	if( strlen( opt_binary_file ) > 0 ) {
		FILE *fp = fopen( opt_binary_file, "rb" );
		if( fp ) {
			if( _save( &m, fp ) )
				errx( -1, "failed saving binary representation of %s", FNAME ? FNAME : "stdin" );
			fclose( fp );
		}
	}
#endif

	// Echo the matrix as text

	if( opt_echo ) {
		if( opt_expect_rownames ) {
			mtm_resolve_rownames( &m, (signed long)m.row_names );
		}
		_dump( &m, opt_label_format, opt_float_format, stdout );
#ifdef HAVE_MD5
		if( opt_include_md5 ) {
			fprintf( stdout, "# MD5: %s\n", m.md5 );
		}
#endif
	}

	m.destroy( &m );

	return errnum ? -1 : 0;
}

