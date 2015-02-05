
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <err.h>
#include <alloca.h>

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
	size_t allocation, tailsize;

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

#if 0
	/**
	  * The header is not persisted (in RAM, as part of the struct 
	  * mtm_matrix), so we can subtract the size of its section from the
	  * tailsize memory.
	  */
	{
		const int LAST_SECTION
			= (header->flags & MTMHDR_ROW_LABELS_PRESENT) != 0
			? S_ROWMAP
			: S_DESCRIPTOR;
		tailsize
			= header->section[ LAST_SECTION ].offset
			+ header->section[ LAST_SECTION ].size
			- header->section[ S_MATRIX ].offset;// because the header need not persist
	}
#else
	tailsize
		= info.st_size - header->section[ S_MATRIX ].offset;
	allocation
		= page_aligned_ceiling( tailsize );
#endif

	if( posix_memalign( (void**)&pc, RT_PAGE_SIZE, allocation ) == 0 ) {

		const size_t SIZEOF_HEADER_SECTION
			= page_aligned_ceiling( header->section[ S_MATRIX ].offset );

		if( fseek( fp, SIZEOF_HEADER_SECTION, SEEK_SET ) )
			return MTM_E_IO;

		if( fread( pc, sizeof(char), tailsize, fp ) != tailsize ) {
			free( pc );
			return MTM_E_IO;
		}

		matrix->data
			= (mtm_int_t*)pc;
		matrix->prop
			= (struct mtm_descriptor *)(pc
			+ header->section[ S_DESCRIPTOR ].offset
			- SIZEOF_HEADER_SECTION);

		if( header->section[ S_ROWID ].offset > 0 ) {
			matrix->row_names
				= (const char *)(pc
				+ header->section[ S_ROWID ].offset
				- SIZEOF_HEADER_SECTION);
		} else
			matrix->row_names = NULL;

		if( header->section[ S_ROWMAP ].offset > 0 ) {
			matrix->row_map
				= (struct mtm_row_id *)(pc
				+ header->section[ S_ROWMAP ].offset
				- SIZEOF_HEADER_SECTION);
		} else
			matrix->row_map = NULL;

		matrix->destroy = mtm_free_matrix;
		matrix->storage = pc;

	} else {
		warnx( "failed allocating RAM in %s", __func__ );
		return MTM_E_NOMEM;
	}

	if( matrix->row_names != NULL && matrix->row_map != NULL )
		mtm_resolve_rownames( matrix, (signed long)matrix->row_names );

	matrix->lexigraphic_order
		= ((header->flags & MTMHDR_ROW_LABELS_LEXORD) != 0);

	return MTM_OK;
}

