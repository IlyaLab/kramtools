
/**
  * See README.rst for format details.
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <sys/sendfile.h>
#include <sys/types.h>   // for lseek
#include <unistd.h>      // for lseek
#include <errno.h>
#include <err.h>
#include <assert.h>

#ifdef HAVE_MD5
#include "md5/md5.h"
#ifndef MD5_DIGEST_LENGTH
#define MD5_DIGEST_LENGTH (16)
#endif
#endif

#include "syspage.h"
#include "mtheader.h"
#include "mtmatrix.h"
#include "feature.h"
#include "mtsclass.h"
#include "mterror.h"
#include "specialc.h"

extern void mtm_free_matrix( struct mtm_matrix *m );

/***************************************************************************
  * Input parameters
  */

/**
  * These are redefinable by defining environment vars:
  * MTM_COMMENT_FLAG
  * MTM_CHAR_FIELD_SEP_CHAR
  */
static const char *ENVVAR_CHAR_FIELD_SEP
	= "MTM_SEPARATOR_CHAR";
static const char *ENVVAR_COMMENT
	= "MTM_COMMENT_CHAR";

/***************************************************************************
  * Parse result/output parameters
  */


static void _pad_to_pagesize( FILE *fp ) {
	if( RT_PAGE_SIZE == 0 )
		page_aligned_ceiling( 1 ); // just to initialize RT_PAGE_SIZE.
	while( ftell(fp) & RT_PAGE_MASK )
		fputc( 0, fp );
}


/**
  * Parse a text matrix satisfying the format description (elsewhere).
  *
  * Parsing means:
  * 1. the file is converted from text to binary form
  * 2. the format is implicitly validated
  * 3. univariate degeneracies in data rows are detected and characterized.
  *
  * This function can leave its results in one (or both) of two forms:
  * 1. memory resident, described by a struct mtm_matrix.
  * 2. file resident, in the caller-provided FILE pointer.
  *
  * Basic operation:
  *
  * This function basically just assumes the caller wants a filesystem-
  * resident result and builds the persistent binary image of the input
  * matrix...either in the caller provided FILE pointer or a tmp file.
  *
  * If the caller wants a RAM-resident result, the built file is then
  * simply loaded by the mtm_load_matrix API.
  *
  * This approach may incur some unnecessary work, but the slight 
  * inefficiency is preferable to maintaining two independent code paths.
  *
  * To minimize gratuitous copying, the matrix' data is cached directly in
  * the caller-provided FILE pointer at its proper final place, so at least
  * the largest block of bytes will not be unnecessarily copied.
  */
int mtm_parse( FILE *input,
		unsigned int flags,
		const char *missing_data_regex,
		int max_allowed_categories,
		int (*infer_stat_class)(const char *),
		FILE *output_fp,
		struct mtm_matrix *output_m ) {

	const size_t SIZEOF_HEADER_BLOCK
		= page_aligned_ceiling(sizeof(struct mtm_matrix_header));
	const int verbosity
		= MTM_VERBOSITY_MASK & flags;
	const bool EXPECT_ROW_NAMES
		= ( flags & MTM_MATRIX_HAS_ROW_NAMES) != 0;
	const bool EXPECT_HEADER
		= ( flags & MTM_MATRIX_HAS_HEADER )   != 0;
	const bool PRESERVE_ROWNAMES
		= EXPECT_ROW_NAMES && ( (flags & MTM_DISCARD_ROW_NAMES) == 0 );

	int   econd = MTM_OK;
	int    lnum = 0;
	int    fnum = 0;

	char  *line = NULL;
	size_t blen = 0;
	ssize_t llen;

	struct feature f = {
		.length = 0,
		//.buf.num: NULL,
		.expect_row_labels =  EXPECT_ROW_NAMES,
		.missing_data_regex = missing_data_regex,
		.interpret_prefix =   infer_stat_class,
		.max_cardinality =    max_allowed_categories,
		.category_labels =    NULL
	};
	struct mtm_descriptor d;

	/**
	  * If caller wants the binary result stored in a file, <output_fp> should
	  * be non-NULL. In this case the data will be written directly into
	  * the file, bypassing tmp files altogether.
	  */

	FILE *data_fp = output_fp ? output_fp : tmpfile();
	// Temp file caches. MATRIX was allocated above this scope.
	// We always need DESCRIPTOR and may or may not need ROW sections.
	FILE *tmp_section[ S_COUNT ];

	struct mtm_matrix_header hdr;

	memset( &hdr, 0, sizeof(struct mtm_matrix_header) );

	if( data_fp == NULL )
		return MTM_E_IO;

	/**
	  * Position the file pointer to leave room for the header which
	  * we'll "back up" and write as the very last step.
	  */

	if( fseek( data_fp, SIZEOF_HEADER_BLOCK, SEEK_SET ) )
		return MTM_E_IO;

	// Preceding was last early return; must execute cleanup at bottom.

	/**
	  * Temp files to contain unpredictably-sized arrays, each of which
	  * corresponds to one of the data sections in struct mtm_matrix_header. 
	  * NOTE THAT THE ONE INDEXED BY S_DATA IS NOT USED.
	  */

	tmp_section[ S_DESC ] = tmpfile();
	if( tmp_section[ S_DESC ] == NULL )
		goto cleanup_tmpfile;

	if( PRESERVE_ROWNAMES ) {

		tmp_section[ S_ROWID  ] = tmpfile();
		if( tmp_section[ S_ROWID ] == NULL )
			goto cleanup_tmpfile;

		tmp_section[ S_ROWMAP ] = tmpfile();
		if( tmp_section[ S_ROWMAP ] == NULL )
			goto cleanup_tmpfile;
	}

#ifdef HAVE_MD5
	md5_state_t hashstate;
	md5_byte_t checksum[ MD5_DIGEST_LENGTH ];
	md5_init( &hashstate );
#endif

	/**
	  * Possibly redefine CHAR_FIELD_SEP and COMMENT flags.
	  */

	if( getenv(ENVVAR_COMMENT) )
		CHAR_COMMENT = getenv(ENVVAR_COMMENT)[0];
	if( getenv(ENVVAR_CHAR_FIELD_SEP) )
		CHAR_FIELD_SEP = getenv(ENVVAR_CHAR_FIELD_SEP)[0];

	/**
	  * Basic parse strategy is line-oriented:
	  * 1) read one line at a time
	  * 2) parse and process it (AFTER checksumming the unadulterated line)
	  * 3) store results in temporary file "caches"
	  * 4) after file is consumed, reconstitute the caches into RAM
	  */
	while( ( llen = getline( &line, &blen, input ) ) > 0 ) {

		++lnum; // ...at beginning of loop for 1-based line reporting.

		assert( line[llen] == 0 ); // ...getline was well-behaved.

		// Checksum before ANY changes happen to the read bytes...

#ifdef HAVE_MD5
		md5_append( &hashstate, (md5_byte_t*)line, llen );
#endif
		// Now replace newline, if present, with NUL terminator to
		// GUARANTEE NUL termination for all subsequent code.
		// (Last line of file might not have a newline.)

		if( line[llen-1] == CHAR_LINE_TERM ) line[--llen] = 0;

		// Skip empty lines and comments.

		if( llen == 0 || line[0] == CHAR_COMMENT )
			continue;

		// We have a non-empty, non-comment line.
		// If it's the first, we have some setup to do...

		if( f.length == 0 ) {

			// Column count of first non-empty, non-comment line establishes
			// the column count for the rest of the file!

			f.length
				= feature_count_fields( line, CHAR_FIELD_SEP )
				- ( EXPECT_ROW_NAMES ? 1 : 0 );

			if( f.length < 2 /* absolute minimum sensible */ ) {
				econd = MTM_E_FORMAT_MATRIX;
				break;
			}
			if( feature_alloc_encode_state( &f ) ) {
				econd = MTM_E_SYS;
				break;
			}

			if( EXPECT_HEADER /* i.e. 1st non-comment line is header */ )
				continue; // no further processing on this line.
		}

		if( ( econd = feature_encode( line, &f, &d ) ) ) {
			warnx( "%s: aborting parsing at input line %d", __FILE__, lnum );
			break;
		}

		/**
		  * Write encoded data, its descriptor and row label data to caches.
		  */

		if( PRESERVE_ROWNAMES ) {
			const struct mtm_row srn = {
				fnum,
				(const char*)ftell(tmp_section[S_ROWID])
			};
			f.label_length++; // ...to include NUL-terminator in the output.
			if( fwrite( line, sizeof(char), f.label_length, tmp_section[S_ROWID] ) != f.label_length
				||
				fwrite( &srn, sizeof(struct mtm_row), 1, tmp_section[S_ROWMAP] )
					!= 1 ) {
				econd = MTM_E_IO;
				break;
			}
		}

		if( fwrite( &d, sizeof(struct mtm_descriptor), 1, tmp_section[S_DESC] )
				!= 1 ) {
			econd = MTM_E_IO;
			break;
		}
		if( fwrite( f.buf.cat, sizeof(mtm_int_t), f.length, data_fp ) != f.length ) {
			econd = MTM_E_IO;
			break;
		}

		fnum += 1;
	}

	if( line )
		free( line );

	if( econd == MTM_OK ) {

#ifdef HAVE_MD5
		md5_finish( &hashstate, checksum );
		for(int i = 0; i < MD5_DIGEST_LENGTH; i++ )
			sprintf( hdr.md5 + 2*i, "%02x", checksum[i] );
#endif

		/**
		  * Build the result in data_fp. If file output was requested, then
		  * data_fp *is* that file. Otherwise, data_fp is a tmpfile.
		  * If file output was not requested--that is, if caller only wants
		  * a RAM-resident matrix--then this is an unnecessary step, but
		  * extra work here is preferable to maintaining two code paths!
		  */

		memcpy( hdr.sig, MTM_SIGNATURE, sizeof(hdr.sig) );
		hdr.endian      = 0x04030201;
		hdr.version     = 0x01000000;
		hdr.flags       = PRESERVE_ROWNAMES ? MTMHDR_ROW_LABELS_PRESENT : 0;
		// The string table is always saved in the matrix' row order, not lexigraphic.
		hdr.header_size = sizeof(struct mtm_matrix_header);
		hdr.sizeof_cell = sizeof(mtm_int_t);
		hdr.rows        = fnum;
		hdr.columns     = f.length;

		/**
		  * The sizes of each section are the current file offsets of the
		  * corresponding caches (except the S_DATA section which needs
		  * the header block subtracted).
		  */

		hdr.section[ S_DATA ].offset
			= SIZEOF_HEADER_BLOCK;
		hdr.section[ S_DATA ].size
			= ftell( data_fp )
			- hdr.section[ S_DATA ].offset;

		// The other sections are simple (no header prefixes).

		hdr.section[ S_DESC ].size
			= ftell( tmp_section[ S_DESC ] );
		rewind( tmp_section[ S_DESC ] );

		// Row label sections may not even be present...

		if( tmp_section[ S_ROWID ] ) {
			hdr.section[ S_ROWID ].size
				= ftell( tmp_section[ S_ROWID ] );
			rewind( tmp_section[ S_ROWID ] );
		}

		if( tmp_section[ S_ROWMAP ] ) {
			hdr.section[ S_ROWMAP ].size
				= ftell( tmp_section[ S_ROWMAP ] );
			rewind( tmp_section[ S_ROWMAP ] );
		}

		// Copy each non-empty section into the file sequentially
		// and beginning on page boundaries.

		for(int i = S_DESC; i < S_COUNT; i++ ) {
			if( tmp_section[ i ] ) {
				_pad_to_pagesize( data_fp );
				hdr.section[ i ].offset = ftell( data_fp );
				sendfile(
					fileno( data_fp ),
					fileno( tmp_section[i] ),
					NULL,
					hdr.section[i].size );

				// Need to keep the C library FILE* in sync with the
				// preceding manipulation of the lower-level file.
				// This is one way to do it...

				fseek( data_fp, 
					lseek( fileno( data_fp ), 0, SEEK_CUR ),
					SEEK_SET );

				assert( lseek(fileno( data_fp ),0,SEEK_CUR) == ftell(data_fp ) );
			}
		}

		/**
		  * Note that padding happens before each section write above, so the
		  * last section does not (and need not) have trailing padding. Thus,
		  * the file size (SIZEOF_HEADER_BLOCK + hdr.sizeof_rt_image) will not
		  * typically be a multiple of PAGE_SIZE.
		  */

		hdr.sizeof_rt_image 
			= ftell( data_fp ) 
			- SIZEOF_HEADER_BLOCK;

		/**
		  * Finally write the header, now that all its contents are complete.
		  */

		rewind( data_fp );
		fwrite( &hdr, sizeof(struct mtm_matrix_header), 1, data_fp );
		_pad_to_pagesize( data_fp );

		/**
		  * If caller wanted a memory-resident result, just reload it.
		  */

		if( output_m ) {
			rewind( data_fp );
			mtm_load_matrix( data_fp, output_m, NULL );
		}
	}

	feature_free_encode_state( &f );

cleanup_tmpfile:

	if( output_fp == NULL /* implying data_fp is a tmpfile */ ) {
		if( data_fp )
			fclose( data_fp );
	}

	for(int i = 0; i < S_COUNT; i++ ) {
		if( tmp_section[i] != NULL )
			fclose( tmp_section[i] );
	}

	return econd;
}

