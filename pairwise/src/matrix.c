
#include <stdlib.h>
#include "matrix.h"

/**
  * Used only by the bsearch function to look up row names.
  */
static int _cmp_row_names( const void *pvl, const void *pvr ) {
	ROW_T *l = (ROW_T*)pvl;
	ROW_T *r = (ROW_T*)pvr;
	return strcmp( l->string, r->string );
}


/**
 * This, on the other hand is used to re-sort the row table for offset
 * lookup. The table arrives sorted by name (actually created in that order
 * by the preprocessor) to facilitate name lookup. If we're running in batch
 * mode we'll convert offsets to names.
 */
static int _cmp_row_rowoffs( const void *pvl, const void *pvr ) {
	ROW_T *l = (ROW_T*)pvl;
	ROW_T *r = (ROW_T*)pvr;
	return l->rowoff==r->rowoff ? 0 : ( l->rowoff < r->rowoff?-1:+1 );
}



/**
  * The row map, upon creation as well as upon loading from a preprocessed file
  * contains offsets from the base of a hypothetical block of memory containing
  * packed, null-terminated strings. Those offsets need to be converted into
  * actual -valid- pointers before analysis time.
  *
  * Depending on the mode of use of the matrix, the table may also need to be
  * (re)sorted.
  */
void rtm_relocate_rowname( struct matrix *m, const char *base, bool lexographic_order ) {

	for(int i = 0; i < m->rows; i++ )
		m->row_map[i].string += base;

	if( m->lexographic_order != lexographic_order ) {
		qsort( g_rowmap, _nrow, sizeof(ROW_T),
			lexographic_order ? _cmp_row_names : _cmp_row_rowoffs );
		m->lexographic_order = lexographic_order;
	}
}

