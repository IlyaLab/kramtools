
/**
  * Yet Another Matrix struct declaration, tailored of course for this 
  * app's requirements.
  *
  * In particular, the runtime may acquire a matrix from two very different sources:
  * 1. a preprocessed file that is memory-mapped.
  * 2. a text file that is preprocessed just-in-time.
  */

#ifndef _matrix_h_
#define _matrix_h_

struct rowname {
	unsigned    rowoff;
	const char *string;
};
typedef struct rowname rowname_t;
typedef const rowname_t ROWNAME_T;

struct matrix {

	int rows;

	/**
	  * This does NOT include the 32-bit feature qualifier prefix on each
	  * feature row.
	  */
	int data_columns;

	/**
	  * This is actually a mix of unsigned int's and float's, 
	  * both of which are not coincidentally 4-bytes in size.
	  */
	unsigned int *data;

	/**
	  * This is always a free-standing array regardless of how the
	  * matrix was acquired.
	  */
	struct rowname *row_map;

	/**
	  * Indicates that the row_map is in lexographic (string) order
	  * as opposed to the order of row numbers in the input matrix.
	  */
	bool lexographic_order;

	/**
	  * This is only valid if the file was JIT preprocessed.
	  * Otherwise, the rownames reside with the main (memory-mapped)
	  * body of the file (and won't need to be indpendently freed).
	  */
	const char *row_names;

	/**
	  * The clean-up function set by whatever created this matrix.
	  */
	void (*destroy)( struct matrix * );

	/**
	  * A place for the creator of this matrix to stash private info
	  * that the destroyer needs to do its job exhaustively.
	  */
	void *aux;
};

void rtm_relocate_rowname( struct matrix *m, const char *base, bool lexographic_order );

#endif

