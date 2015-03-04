
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <err.h>
#include <alloca.h>
#include <assert.h>

#include "mtmatrix.h"
#include "mtheader.h"
#include "mterror.h"
#include "syspage.h"

/**
  * The parser allocates the matrix as one large memory blob, so
  * this destroyer is particularly simple.
  */
void mtm_free_matrix( struct mtm_matrix *m ) {
	if( m ) {
		if( m->storage )
			free( m->storage );
		m->storage = NULL;
	}
}


int mtm_load_header( FILE *fp, struct mtm_matrix_header *header ) {

	if( header == NULL || fp == NULL )
		return MTM_E_NULLPTR;

	if( fread( header, sizeof(struct mtm_matrix_header), 1, fp ) != 1 )
		return MTM_E_IO;

	if( strcmp( header->sig, MTM_SIGNATURE ) )
		return MTM_E_BADSIG;

	return MTM_OK;
}


int mtm_load_matrix( FILE *fp, struct mtm_matrix *matrix, struct mtm_matrix_header *header ) {

	char *pc = NULL;
	struct stat info;
	size_t allocation;

	if( fstat( fileno( fp ), &info ) )
		return MTM_E_IO;

	if( matrix == NULL )
		return MTM_E_NULLPTR;

	if( header == NULL )
		header = alloca( sizeof(struct mtm_matrix_header) );

	if( fread( header, sizeof(struct mtm_matrix_header), 1, fp ) != 1 )
		return MTM_E_IO;

	matrix->rows    = header->rows;
	matrix->columns = header->columns;
	matrix->size    = header->sizeof_rt_image;

	assert( matrix->size 
			= info.st_size 
			- header->section[ S_DATA ].offset );

	allocation // ...must be PAGE_SIZE multiple, so will be a bit more.
		= page_aligned_ceiling( matrix->size );

	if( posix_memalign( (void**)&pc, RT_PAGE_SIZE, allocation ) == 0 ) {

		const size_t SIZEOF_HEADER_BLOCK
			= page_aligned_ceiling(sizeof(struct mtm_matrix_header));

		memset( pc, 0, allocation );

		if( fseek( fp, SIZEOF_HEADER_BLOCK, SEEK_SET ) )
			return MTM_E_IO;

		if( fread( pc, sizeof(char), matrix->size, fp ) != matrix->size ) {
			free( pc );
			return MTM_E_IO;
		}

		matrix->data
			= (mtm_int_t*)pc;
		matrix->desc
			= (struct mtm_descriptor *)(pc
			+ header->section[ S_DESC ].offset
			- SIZEOF_HEADER_BLOCK);

		if( header->section[ S_ROWID ].offset > 0 ) {
			matrix->row_id
				= (const char *)(pc
				+ header->section[ S_ROWID ].offset
				- SIZEOF_HEADER_BLOCK);
		} else
			matrix->row_id = NULL;

		if( header->section[ S_ROWMAP ].offset > 0 ) {
			matrix->row_map
				= (struct mtm_row *)(pc
				+ header->section[ S_ROWMAP ].offset
				- SIZEOF_HEADER_BLOCK);
		} else
			matrix->row_map = NULL;

		matrix->destroy = mtm_free_matrix;
		matrix->storage = pc;

	} else {
		warnx( "failed allocating RAM in %s", __func__ );
		return MTM_E_NOMEM;
	}

	if( matrix->row_id != NULL && matrix->row_map != NULL )
		mtm_resolve_rownames( matrix, (signed long)matrix->row_id );

	matrix->lexigraphic_order
		= ((header->flags & MTMHDR_ROW_LABELS_LEXORD) != 0);

	return MTM_OK;
}

