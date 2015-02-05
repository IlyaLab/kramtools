
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <err.h>

#include "syspage.h"
#include "mtmatrix.h"
#include "mtheader.h"
#include "mterror.h"


int main( int argc, char *argv[] ) {

	struct mtm_matrix_header hdr;
	FILE *fp = NULL;
	const char *fname = NULL;
	int econd;

	if( argc < 2 ) {
		printf( "%s <preprocessed matrix>\n", argv[0] );
		exit(-1);
	} else
		fname = argv[1];

	fp = fopen( fname, "r" );
	if( NULL == fp )
		err( -1, "opening %s", fname );

	if( (econd = mtm_load_header( fp, &hdr )) == MTM_OK ) {

		int lnum = 1;
		char *line = NULL;
		size_t blen = 0;
		ssize_t llen = 0;

		const size_t SIZEOF_STRINGS
			= hdr.section[ S_ROWID ].size;
		const size_t SIZEOF_MAP
			= hdr.section[ S_ROWMAP ].size;
		const char *strings
			= calloc( SIZEOF_STRINGS, sizeof(char) );
		struct mtm_row_id *map
			= calloc( hdr.rows, sizeof(struct mtm_row_id) );

		if( strings == NULL || map == NULL ) {
			err( -1, "failed allocating one or both of %ld and %d bytes",
				SIZEOF_STRINGS, SIZEOF_MAP );
		}

		/**
		  * Load the strings...
		  */

		if( fseek( fp, hdr.section[ S_ROWID ].offset, SEEK_SET ) )
			err( -1, "failed seeking to string table" );
		if( fread( (void *)strings, sizeof(char), SIZEOF_STRINGS, fp ) != SIZEOF_STRINGS )
			err( -1, "failed loading string table" );

		/**
		  * ...and the map.
		  */

		if( fseek( fp, hdr.section[ S_ROWMAP ].offset, SEEK_SET ) )
			err( -1, "failed seeking to row map" );
		if( fread( (void *)map, sizeof(struct mtm_row_id ), hdr.rows, fp ) != hdr.rows )
			err( -1, "failed loading row map" );

		fclose( fp );
		fp = NULL;

		/**
		  * Update the map to point into the string table.
		  */

		for(int i = 0; i < hdr.rows; i++ )
			map[i].string += (long)strings;

		while( ( llen = getline( &line, &blen, stdin ) ) > 0 ) {
			const int ROW
				= atoi( line );
			if( 0 <= ROW && ROW < hdr.rows )
				puts( map[ROW].string );
			else {
				fprintf( stderr, "%d on line %d is invalid row index in matrix with %d rows\n",
					ROW, lnum, hdr.rows );
			}
			lnum ++;
		}
		if( line )
			free( line );
fail:
		if( strings )
			free( strings );
		if( map )
			free( map );

	} else
	if( econd == MTM_E_BADSIG )
		errx( -1, "%s has wrong signature.\n"
			"Are you sure this is a preprocessed matrix?",
			fname );
	else
		err( -1, "failed loading header from %s", fname );

	if( fp )
		fclose( fp );

	return EXIT_SUCCESS;
}

