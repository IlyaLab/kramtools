
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "mtmatrix.h"
#include "mterror.h"

/**
  * Used only by the bsearch function to look up row names.
  */
static int _cmp_key_and_row_id( const void *pvl, const void *pvr ) {
	const char *l
		= (const char *)pvl;
	struct mtm_row *r
		= (struct mtm_row*)pvr;
	return strcmp( l, r->string );
}


/**
  * Used only by the qsort to sort by row names.
  */
static int _cmp_row_id( const void *pvl, const void *pvr ) {
	struct mtm_row *l = (struct mtm_row*)pvl;
	struct mtm_row *r = (struct mtm_row*)pvr;
	return strcmp( l->string, r->string );
}


/**
 * This, on the other hand is used to re-sort the row table for offset
 * lookup. The table arrives sorted by name (actually created in that order
 * by the preprocessor) to facilitate name lookup. If we're running in batch
 * mode we'll convert offsets to names.
 */
static int _cmp_mtm_row_offsets( const void *pvl, const void *pvr ) {
	struct mtm_row *l = (struct mtm_row*)pvl;
	struct mtm_row *r = (struct mtm_row*)pvr;
	return l->offset==r->offset ? 0 : ( l->offset < r->offset?-1:+1 );
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
void mtm_resolve_rownames( struct mtm_matrix *m, signed long base ) {
	if( m->row_map ) {
		for(int i = 0; i < m->rows; i++ )
			m->row_map[i].string += base;
	}
}


int mtm_resort_rowmap( struct mtm_matrix *m, bool lexigraphic_order ) {
	if( m->row_map ) {
		if( m->lexigraphic_order != lexigraphic_order ) {
			qsort( m->row_map, m->rows, sizeof(struct mtm_row),
				lexigraphic_order ? _cmp_row_id : _cmp_mtm_row_offsets );
			m->lexigraphic_order = lexigraphic_order;
		}
		return 0;
	}
	return MTM_E_NO_ROW_LABELS;
}


#ifdef _DEBUG

/**
  * Following are sheer paranoia...
  */

static bool _sorted_lexigraphic( MTM_ROW_ID_T *map, int n ) {
	for(int i = 1; i < n; i++ ) {
		if( strcmp( map[i-1].string, map[i].string ) > 0 )
			return false;
	}
	return true;
}

static bool _sorted_by_row_offset( MTM_ROW_ID_T *map, int n ) {
	for(int i = 1; i < n; i++ ) {
		if( map[i-1].offset > map[i].offset )
			return false;
	}
	return true;
}

#endif

/**
  * mtm_fetch_by_name and mtm_fetch_by_offset are essentially convenience 
  * functions. They gather into the struct mtm_feature all the (scattered)
  * elements of one feature, relieving the libmtm client of dealing directly
  * with the struct mtm_matrix (even though it's not opaque).
  */
int  mtm_fetch_by_name( struct mtm_matrix *m, struct mtm_feature *f ) {

	struct mtm_row *prid;

	assert( m->lexigraphic_order );
	assert( _sorted_lexigraphic( m->row_map, m->rows ) );

	if( m->row_map == NULL )
		return MTM_E_NO_ROW_LABELS;

	if( ( m == NULL ) 
			|| ( f == NULL ) 
			|| ( f->name == NULL ) )
		return MTM_E_NULLPTR;

	prid = bsearch( f->name, m->row_map, m->rows, sizeof(struct mtm_row), _cmp_key_and_row_id );
	if( prid ) {
		const int ROW
			= prid->offset;
		f->offset = ROW;
		f->desc   = m->desc[ ROW ];
		f->data   = m->data + ROW*m->columns;
		return 0;
	}
	return MTM_E_NO_SUCH_FEATURE;
}


int  mtm_fetch_by_offset( struct mtm_matrix *m, struct mtm_feature *f ) {

	const int ROW
		= f->offset;

	if( m == NULL || f == NULL )
		return MTM_E_NULLPTR;

	assert( ! m->lexigraphic_order );
	assert( m->row_map == NULL || _sorted_by_row_offset( m->row_map, m->rows ) );

	if( 0 <= ROW && ROW < m->rows ) {
		f->name = m->row_map ? m->row_map[ ROW ].string : NULL;
		f->desc = m->desc[ ROW ];
		f->data = m->data + ROW*m->columns;
		return 0;
	}
	return MTM_E_NO_SUCH_FEATURE;
}

