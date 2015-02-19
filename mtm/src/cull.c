
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

#include "md5/md5.h"
#ifndef MD5_DIGEST_LENGTH
#define MD5_DIGEST_LENGTH (16)
#endif


struct feature_hash {
	md5_byte_t digest[ MD5_DIGEST_LENGTH ];
	unsigned offset;
};

static int _cmp_feature_hash( const void *pvl, const void *pvr ) {

	const struct feature_hash *l
		= (const struct feature_hash *)pvl;
	const struct feature_hash *r
		= (const struct feature_hash *)pvr;

	return memcmp( l->digest, r->digest, MD5_DIGEST_LENGTH );
}


struct block {
	unsigned start,end;
};


static int _cmp_block_r( const void *pvl, const void *pvr ) {
	const struct block *l
		= (const struct block *)pvl;
	const struct block *r
		= (const struct block *)pvr;
	return -( (int)(l->end - l->start) - (int)(r->end - r->start) );
}


static int _load_feature_names( FILE *fp, struct mtm_matrix_header *hdr,
	struct mtm_row_id **map,
	const char **strings ) {

	const size_t SIZEOF_STRINGS
		= hdr->section[ S_ROWID ].size;
	const size_t SIZEOF_MAP
		= hdr->section[ S_ROWMAP ].size;
	*strings
		= calloc( SIZEOF_STRINGS, sizeof(char) );
	*map
		= calloc( hdr->rows, sizeof(struct mtm_row_id) );

	if( *strings == NULL || *map == NULL ) {
		err( -1, "failed allocating one or both of %ld and %d bytes",
			SIZEOF_STRINGS, SIZEOF_MAP );
	}

	/**
	  * Load the strings...
	  */

	if( fseek( fp, hdr->section[ S_ROWID ].offset, SEEK_SET ) )
		err( -1, "failed seeking to string table" );
	if( fread( (void*)(*strings), sizeof(char), SIZEOF_STRINGS, fp ) != SIZEOF_STRINGS )
		err( -1, "failed loading string table" );

	/**
	  * ...and the map.
	  */

	if( fseek( fp, hdr->section[ S_ROWMAP ].offset, SEEK_SET ) )
		err( -1, "failed seeking to row map" );
	if( fread( *map, sizeof(struct mtm_row_id ), hdr->rows, fp ) != hdr->rows )
		err( -1, "failed loading row map" );

	/**
	  * Update the map to point into the string table.
	  */

	for(int i = 0; i < hdr->rows; i++ )
		(*map)[i].string += (long)(*strings);

	return 0;
}


static int _report( FILE *blocks, FILE *matrix, const struct feature_hash *hashes, const struct mtm_matrix_header *hdr ) {

	const size_t BLOCKS
		= ftell( blocks ) / sizeof(struct block);
	struct block *b
		= calloc( BLOCKS, sizeof(struct block) );
	struct mtm_row_id *map = NULL;
	const char *strings    = NULL;

	if( b == NULL )
		err( -1, "failed allocating array of blocks" );

	/**
	  * Load and sort the blocks by size
	  */
	rewind( blocks );

	if( fread( b, sizeof(struct block), BLOCKS, hashes ) == BLOCKS ) {
		qsort( h, BLOCKS, sizeof(struct block), _cmp_block );
	} else
		err( -1, "failed rereading blocks file" );

	if( _load_feature_names( fp, hdr, &map, &strings ) == 0 )

		if( strings )
			free( (void*)strings );
		if( map )
			free( (void*)map );
	} else
		err( -1, "failed loading feature names and map" );
}


/**
  * Reloads the hashes from the tmpfile, sorts them then identifies
  * contiguous ranges of identical hashes (indicating identical features).
  * Writes these to another tmp file.
  */
static int _analyze( FILE *hashes, FILE *fp, const struct mtm_matrix_header *hdr ) {

	const size_t HASHED_FEATURES
		= ftell( hashes ) / sizeof(struct feature_hash);
	struct feature_hash *h
		= calloc( HASHED_FEATURES, sizeof(struct feature_hash) );
	FILE *blocks = NULL;

	unsigned n = 1, start = 0;
	if( h == NULL )
		return -1;

	/**
	  * Load and sort (and maybe rewrite to the same tmp file) the feature hashes.
	  */

	rewind( hashes );

	if( fread( h, sizeof(struct feature_hash), HASHED_FEATURES, hashes ) == HASHED_FEATURES ) {
		qsort( h, HASHED_FEATURES, sizeof(struct feature_hash), _cmp_feature_hash );
#if 0
		ftruncate( fileno(hashes), 0L );
		if( fwrite( h, sizeof(struct feature_hash), HASHED_FEATURES, hashes ) != HASHED_FEATURES )
			err( -1, "failed rewriting sorted feature hashes to tmp file" );
		free( h ); h = NULL;
#endif
	} else
		err( -1, "failed rereading hashed features file" );

	/**
	  * Now identify contiguous blocks of identical digests indicating
	  * (almost certainly) redundant features.
	  */

	blocks = tmpfile();
	if( blocks ) {
		n = 1, start = 0;
		for(int i = 1; i < HASHED_FEATURES; i++ ) {
			if( memcmp( h[start].digest, h[i].digest, MD5_DIGEST_LENGTH ) ) {
				if( n > 1 ) {
					struct block b = { 
						.start = start,
						.end   = i
					};
					if( fwrite( &b, sizeof(struct block), 1, blocks ) != 1 )
						err( -1, "failed writing a block" );
				}
				start = i;
				n     = 1;
			}	
		}

		_report( blocks, fp, hdr );

		fclose( blocks );

	} else
		err( -1, "failed creating tmp file for blocks" );
	/**
	  * Free everything.
	  */

	if( h )
		free( (void*)h );
}


int main( int argc, char *argv[] ) {

	struct mtm_matrix_header hdr;
	FILE *fpd = NULL, *fpm = NULL;
	const char *fname = NULL;
	int econd;

	if( argc < 2 ) {
		printf( "%s <preprocessed matrix>\n", argv[0] );
		exit(-1);
	} else
		fname = argv[1];

	fpd = fopen( fname, "r" );
	fpm = fopen( fname, "r" );

	if( NULL == fpd )
		err( -1, "opening %s", fname );

	if( (econd = mtm_load_header( fpd, &hdr )) == MTM_OK ) {

		FILE *hashfile = tmpfile();
		const size_t SIZEOF_FEATURE
			= hdr.columns*sizeof(mtm_int_t);
		struct feature_hash f;
		struct mtm_descriptor d;
		void *r
			= calloc( hdr.columns, sizeof(mtm_int_t) );

		/**
		  * Seek to the start of each section.
		  */

		if( fseek( fpd, hdr.section[ S_DESCRIPTOR ].offset, SEEK_SET ) )
			err( -1, "failed seeking to descriptor table" );

		if( fseek( fpm, hdr.section[ S_MATRIX ].offset, SEEK_SET ) )
			err( -1, "failed seeking to descriptor table" );

		/**
		  * Now iterate across the rows of the matrix hashing each
		  * categorical or boolean and save the index and hash.
		  */

		//if( opt_verbosity > 1 )
			puts( "reading matrix..." );

		for(f.offset = 0; f.offset < hdr.rows; f.offset++ ) {

			if( fread( &d, sizeof(struct mtm_descriptor), 1, fpd ) != 1 )
				err( -1, "failed reading descriptor %d", f.offset );

			if( fread( &r, sizeof(mtm_int_t), hdr.columns, fpm ) != hdr.columns )
				err( -1, "failed reading row %d", f.offset );

			if( d.categories > 0 && d.constant == 0 ) {
				md5_state_t h;
				md5_init( &h );
				md5_append( &h, (md5_byte_t*)r, SIZEOF_FEATURE );
				md5_finish( &h, f.digest );
				if( fwrite( &f, sizeof(f), 1, hashfile ) != 1 )
					err( -1, "failed caching hash of feature %d", f.offset );
			}
		}

		//if( opt_verbosity > 1 )
			puts( "analyzing..." );

		free( r ); r = NULL;

		_analyze( hashfile, fpd, &hdr );

	} else
	if( econd == MTM_E_BADSIG )
		errx( -1, "%s has wrong signature.\n"
			"Are you sure this is a preprocessed matrix?",
			fname );
	else
		err( -1, "failed loading header from %s", fname );

	if( fpm )
		fclose( fpm );
	if( fpd )
		fclose( fpd );

	return EXIT_SUCCESS;
}

