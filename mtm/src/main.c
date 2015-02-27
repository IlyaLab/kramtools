
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>    // isnan
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <err.h>

#include "syspage.h"
#include "mtmatrix.h"
#include "mtheader.h"
#include "mterror.h"

extern int mtm_sclass_by_prefix( const char *token );

static void _echo_header( struct mtm_matrix_header *h, FILE *fp ) {
	fprintf( fp, "%s\n", h->sig );
	fprintf( fp,
		"     endian: %08x\n"
		"    version: %08x\n"
		"      flags: %08x\n"
		"header_size: %d\n"
		" datum_size: %d\n"
		"       rows: %d\n"
		"    columns: %d\n",
		h->endian,
		h->version,
		h->flags,
		h->header_size,
		h->datum_size,
		h->rows,
		h->columns );
		
	for(int i = 0; i < S_COUNT; i++ ) {
		const struct section_descriptor *s = h->section + i;
		fprintf( fp, "%d: %016lx bytes @ %016lx\n",
			i,
			s->size,
			s->offset );
	}
}


static void _echo_matrix( struct mtm_matrix *m, 
		const char *label_format,
		const char *float_format, 
		FILE *fp ) {

	const char *MISSING_DATA = "NA";

	MTM_INT_T *pdata = m->data;

	// Start echoing.
	// We never echo the header; it's not preserved in parsing.

	for(int r = 0; r < m->rows; r++ ) {

		const int MISSING
			= m->desc[r].missing;
		fprintf( fp, "%s%s%c:%c:%d:%d",
			m->row_map ? m->row_map[r].string : "",
			m->row_map ? "\t"                 : "",
			MISSING < m->columns ? "FI"[ m->desc[r].integral ] : '?',
			m->desc[r].constant ? '!' : '-',
			m->desc[r].categories,
			MISSING );

		if( m->desc[r].integral ) {

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
static bool opt_echo_matrix     = false;
static bool opt_echo_header     = false;
// This is the minimum precision that random.rb generates which
// is the MOST we can assume is available...
static const char *opt_float_format = "%.1e";
static const char *opt_label_format = "%s\t%d:%d:%d:%d";
static int opt_verbosity = 0;

static const char *USAGE = 
"%s [ options ] [ <input file> ] [ <output file> ]\n"
"Input options:\n"
"   --nolabels | -r   do NOT expect input to have row names \n"
"   --noheader | -h   do NOT expect input to have a header \n"
"   --missing  | -m   set the \"missing data\" regex [ \"%s\" ]\n"
"   --maxcats  | -k   set the maximum number of categories allowed \n"
"                     categorical variables [%d]\n"
"   --infer    | -i   infer statistical class from syntax\n"
#if 0
"                     Uses prefix convention describe below by default.\n"
"                     \"C:\" categorical\n"
"                     \"B:\" boolean\n"
"                     \"N:\" numeric/floating-point data\n"
#endif
#ifdef HAVE_MD5
"  --checksum  | -c   include the (MD5) checksum of the input[%s]\n"
#endif
"Testing options (for validating a preprocessed matrix):\n"
"  --echo      | -E   \"echo\" a preprocessed matrix as text [%s]\n"
"  --header    | -H   \"echo\" just the header of a preprocessed matrix [%s]\n"
"  --float     | -F   a C printf format string for floating point\n"
"                     display [\"%s\"]\n"
"  --label     | -L   a C printf format string for display of the\n"
"                     row name and descriptor (see Notes). [\"%s\"]\n"
"  --verbosity | -v verbosity [%d]\n"
"Notes:\n"
"1. The output filename is only optional with \"echoing\" an already-\n"
"   preprocessed matrix (for testing). Otherwise, an output filename\n"
"   is required.\n"
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
		opt_echo_matrix ? _T : _F,
		opt_echo_header ? _T : _F,
		opt_float_format,
		opt_label_format,
		opt_verbosity );
}


int main( int argc, char *argv[] ) {

	const char *STDIN  = "stdin";
	const char *STDOUT = "stdout";
	int errnum;

	const char *fname_i = STDIN;
	FILE          *fp_i = stdin;
	const char *fname_o = STDOUT;
	FILE          *fp_o = stdout;

	opt_missing_marker = mtm_default_NA_regex;

	if( argc < 2 ) {
		_print_usage( argv[0], stdout );
		exit(0);
	}

	do {
		static const char *CHAR_OPTIONS 
#ifdef HAVE_MD5
			= "rhm:k:icEHF:L:v:?";
#else
			= "rhm:k:iEHF:L:v:?";
#endif
		static struct option LONG_OPTIONS[] = {

			{"nolabels",   0,0,'r'},
			{"noheader",   0,0,'h'},
			{"missing",    1,0,'m'},
			{"maxcats",    1,0,'k'},
			{"infer",      0,0,'i'},
#ifdef HAVE_MD5
			{"checksum",   0,0,'c'},
#endif
			{"echo",       0,0,'E'},
			{"header",     0,0,'H'},
			{"float",      1,0,'F'},
			{"label",      1,0,'L'},

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
		case 'i': opt_interpret_type  = NULL;         break;
#ifdef HAVE_MD5
		case 'c': opt_include_md5     = true;         break;
#endif
		case 'E': opt_echo_matrix     = true;         break;
		case 'H': opt_echo_header     = true;         break;
		case 'F': opt_float_format = optarg;          break;
		case 'L': opt_label_format = optarg;          break;

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

	/**
	  * Infer input/output from the first 0-2 positional arguments on the
	  * assumption that
	  * 1) if both are present they're ordered as <input> <output>
	  * 2) <input> must exist (and be readable)
	  * 3) <output> must not exist
	  */

	switch( argc - optind ) {
	case 0: // input: stdin, output: stdout
		break;
	case 1: // i/o depends on whether argv[optind] names an existing file
		if( access( argv[ optind ], F_OK ) == 0 ) {
			fname_i = argv[optind++ ];
			if( access( fname_i, R_OK ) )
				errx( -1,
					"assuming %s is input because it exists, but it is unreadable.",
					fname_i );
			else
			if( opt_verbosity > 0 )
				warnx( "using \"%s\" as input", fname_i );
		} else { // file doesn't exist so assume it's an output filename
			fname_o = argv[ optind++ ];
			warnx( "using \"%s\" as output", fname_o );
		}
		break;

	case 2:
	default: // More than two args ignore all but next two.
		assert( argc - optind >= 2 );
		fname_i = argv[ optind++ ];
		if( access( fname_i, F_OK ) == 0 ) {
			if( access( fname_i, R_OK ) )
				errx( -1,
					"assuming %s is input because it exists, but it is unreadable.",
					fname_i );
		} else
			errx( -1,
				"command line position of \"%s\" implies it's your input, but it doesn't exist.",
				fname_i );
		fname_o = argv[ optind++ ];
		if( access( fname_o, F_OK ) == 0 )
			errx( -1,
				"command line position of \"%s\" implies it's your output, but it exists.\nWon't overwrite.",
				fname_i );
	}

	/**
	  * Open any actual files inferred above.
	  */

	if( strcmp( fname_o, STDOUT ) ) {
		fp_o = fopen( fname_o, "wb" );
	} else
	if( ! ( opt_echo_matrix || opt_echo_header ) )
		errx( -1, "an output filename is not optional when preprocessing" );

	if( strcmp( fname_i, STDIN ) )
		fp_i = fopen( fname_i, "r" );
	else
	if( opt_verbosity > 0 ) // ...so that naive user doesn't assume a hung program.
		warnx( "expecting input on %s\n", STDIN );

	/**
	  * Head must be handled as special case because matrix may be
	  * too large to be loadable in general.
	  */
	if( opt_echo_header ) {
		struct mtm_matrix_header hdr;
		memset( &hdr, 0, sizeof(hdr) );
		if( fread( &hdr, sizeof(struct mtm_matrix_header), 1, fp_i ) == 1 ) {
			_echo_header( &hdr, fp_o );
		} else
			warnx( "failed loading header" );
		rewind( fp_i );
	}

	if( opt_echo_matrix ) {

		struct mtm_matrix mat;
		struct mtm_matrix_header hdr;
		memset( &mat, 0, sizeof(mat) );
		memset( &hdr, 0, sizeof(hdr) );

		if( MTM_OK == mtm_load_matrix( fp_i, &mat, &hdr ) ) {

			_echo_matrix( &mat, opt_label_format, opt_float_format, fp_o );
#ifdef HAVE_MD5
			if( opt_include_md5 ) {
				fprintf( stdout, "# MD5: %s\n", m.md5 );
			}
#endif
			mat.destroy( &mat );
		} else
			warnx( "failed loading %s", fname_i );

	} else {

		const unsigned int FLAGS
			= ( MTM_VERBOSITY_MASK & opt_verbosity)
			| ( opt_expect_rownames ? MTM_MATRIX_HAS_ROW_NAMES : 0 )
			| ( opt_expect_header   ? MTM_MATRIX_HAS_HEADER    : 0 );
		errnum = mtm_parse( fp_i, 
			FLAGS, 
			opt_missing_marker, 
			opt_max_categories, 
			opt_interpret_type, // may be NULL
			fp_o,
			NULL );
	}

	if( fp_o ) fclose( fp_o );
	if( fp_i ) fclose( fp_i );

	return errnum ? -1 : 0;
}

