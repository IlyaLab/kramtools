
#ifndef _mtmatrix_h_
#define _mtmatrix_h_

/**
  * These structs collectively contain an exploded Mixed-Type Matrix.
  * This is a binary matrix containing data elements of uniform -size-, but
  * multiple statistical feature classes:
  *
  * 1. boolean
  * 1. categorical
  * 2. ordinal
  * 3. continuous
  *
  * Rows in this matrix correspond to features; their scalar components to
  * individual samples.
  *
  * By "exploded" I mean it is stored, after parsing, in multiple parts:
  *  1. the data matrix (converted to packed binary form)
  *  2. a vector of row (feature) descriptors (see below)
  * [3. a vector of ( row name, index ) pairs, i.e. the "rowmap"]
  * [4. a packed string table into which the rowmap's string pointers point.]
  *
  * The first 2 items are always present; the row map is optional.
  * The rowmap allows lookup of rows by name or by index depending on
  * how it is sorted.
  */

#ifdef __cplusplus
extern "C" {
#endif

/**
  * Forward declarations
  */
struct mtm_matrix_header;

#ifdef HAVE_LARGE_FORM
typedef unsigned long   mtm_int_t;
typedef const mtm_int_t MTM_INT_T;
typedef double          mtm_fp_t;
typedef const mtm_fp_t  MTM_FP_T;
typedef MTM_INT_T      *MTM_ROW_PTR;
#else
typedef unsigned int    mtm_int_t;
typedef const mtm_int_t MTM_INT_T;
typedef float           mtm_fp_t;
typedef const mtm_fp_t  MTM_FP_T;
typedef MTM_INT_T      *MTM_ROW_PTR;
#endif

/**
  *                          !!! WARNING !!!
  *
  * The integrity of the parsed result hinges on some assumptions about the
  * actual data that this library is likely to encounter and process.
  * In particular, I assume that it is "safe" to encode missing integral
  * data ("NAs") by the same *bit pattern* as (float)NaN (#defined below).
  *
  * For integer-encoded *categorical* data this is perfectly safe since
  * statistical considerations preclude categorical data from involving
  * more than (typically) 10's of categories. In other words, 0x7FC00000
  * is never likely to be needed to represent a category ("factor level").
  *
  * For ordinal data it's still safe, but maybe slightly less so.
  * Ordinal *should* rarely involve 2,143,289,344 levels (the decimal value
  * of 0x7FC00000), but I can't say that it can't.
  *
  * As a result of these considerations, I use the same *bit pattern* (NaN)
  * to represent missing data for all data classes (numeric, ordinal, and
  * categorical) and both data types (floating-point and integral).
  *
  *                          !!! WARNING !!!
  */

#ifdef HAVE_LARGE_FORM
#error "NAN_AS_UINT undefined for 8-byte data sizes"
#else
#define NAN_AS_UINT (0x7FC00000)
#endif

/**
  * This struct stores information about a feature's "quality,"
  */
struct mtm_descriptor {

	unsigned int unused:29;

	/**
	  * For both integral and floating-point features this is the primary
	  * degeneracy flag. It indicates either:
	  * 1. All present fields' values were equal.
	  *    That is, only one unique value was present in the feature.
	  * 2. All fields' values were missing ("NA").
	  *    That is, the feauture is entirely empty of data.
	  * 3. Exactly one value was non-missing (a special case of #1).
	  *
	  * The possibilities are distinguished by looking at <missing>.
	  * If missing == column count, it's case #2.
	  *
	  * In any case, this bit being set indicates the feature is probably
	  * unuseable.
	  */
	unsigned int constant:1;

	/**
	  * This flag indicates the encoding of the data in the feature: either
	  * unsigned int or float. (Does NOT say anything about the
	  * interpretation; it could be boolean, categorical or ordinal.)
	  */
	unsigned int integral:1;

	/**
	  * Says that the data of this feature is categorical (which includes
	  * boolean). Cardinality distinguishes categorical from boolean.
	  */
	unsigned int categorical:1;

	/**
	  * This is only meaningful if integral is set above and it is bounded
	  * by the maximum allowed cardinality+1 of categorical data, so it will
	  * not be accurate if actual data exceeds that.
	  */
	unsigned short cardinality;

	/**
	  * For all feature types, a count of the missing values.
	  *
	  *                          !!! WARNING !!!
	  *
	  * This implies no feature may have more than 65536 missing values.
	  * ...fine for current work, but this limit needs to be highlighted
	  * and asserted far and wide.
	  *
	  *                          !!! WARNING !!!
	  * 
	  */
	unsigned short missing;
};
typedef struct mtm_descriptor mtm_descriptor_t;
typedef const mtm_descriptor_t MTM_DESCRIPTOR_T;
typedef mtm_descriptor_t *mtm_descriptor_ptr;

#define MTM_MAX_MISSING_VALUES (65535)

////////////////////////////////////////////////////////////////////////////

struct mtm_row {
	/**
	  * This is the offset -in the matrix itself- (which is NEVER reordered).
	  */
	unsigned    offset;
	const char *string;
};
typedef struct mtm_row mtm_row_t;
typedef const mtm_row_t MTM_ROW_ID_T;

////////////////////////////////////////////////////////////////////////////

/**
  * This struct is deliberately NOT opaque so that libmtm user can (e.g. in
  * performance critical contexts) iterate the contents of its various 
  * arrays directly (rather than using the arithmetic-heavy mtm_fetch
  * convenience functions).
  */
struct mtm_matrix {

	int rows;
	int columns;

	/**
	  * Assuming the following four sections are contained in a single
	  * block of memory (which need not necessarily be the case), this
	  * is the (minimum) size in bytes of that memory block. (As with
	  * malloc, there may be more allocated than we actually "own.")
	  */
	size_t size;

	/**
	  * This is actually a mix of unsigned int's and float's,
	  * both of which are not coincidentally 4-bytes in size.
	  */
	mtm_int_t *data;

	struct mtm_descriptor *desc;

	/**
	  * This is only valid if the file was JIT preprocessed.
	  * Otherwise, the rownames reside with the main (memory-mapped)
	  * body of the file (and won't need to be indpendently freed).
	  */
	const char *row_id;

	/**
	  * This either maps 
	  * 1. row offsets => row names or 
	  * 2. row names => row offsets
	  * ...depending on how it is sorted (lexigraphic_order). 
	  */
	struct mtm_row *row_map;

	/**
	  * Indicates that the row_map is in lexigraphic (string) order
	  * as opposed to the order of row numbers in the input matrix.
	  */
	bool lexigraphic_order;

	/**
	  * The clean-up function set by whatever filled in this matrix.
	  * This function is responsible for destroying the -contents-
	  * of this structure, NOT freeing the structure itself (which
	  * may be stack-resident!).
	  */
	void (*destroy)( struct mtm_matrix * );

	/**
	  * A place for the creator of this matrix to stash private info
	  * that destroy() needs to do its job exhaustively.
	  * This pair of fields allow for a variety of possible allocation
	  * methods to be used to create a matrix. e.g. A file may be memory
	  * mapped in which case storage may point to a small heap-allocated
	  * struct containing a file handle, offset, and length info.
	  */
	void *storage;
};

void mtm_resolve_rownames( struct mtm_matrix *m, signed long base );

#define MTM_RESORT_LEXIGRAPHIC (true)
#define MTM_RESORT_BYROWOFFSET (false)

int  mtm_resort_rowmap( struct mtm_matrix *m, bool lexigraphic_order );

/**
  * First nibble of flags is used for verbosity.
  */
#define MTM_VERBOSITY_MASK       (0x0000000F)
#define MTM_MATRIX_HAS_HEADER    (0x00000010)
#define MTM_MATRIX_HAS_ROW_NAMES (0x00000020)
#define MTM_DISCARD_ROW_NAMES    (0x00000040)

typedef int (*MTM_ROW_LABEL_INTERPRETER)(const char * );

int mtm_parse( FILE *input,
		unsigned int flags /* first nibble of flags for verbosity */,
		const char *missing_data_regex,
		int max_allowed_categories,
		MTM_ROW_LABEL_INTERPRETER,
		FILE *fout, // may be null
		struct mtm_matrix *m);

/**
  * The immediate motivation for this library is pairwise analysis of
  * features, so the following are convenience APIs for this use case.
  */
struct mtm_feature {
	/**
	  * This is the offset -in the matrix itself- (which is NEVER reordered).
	  * The array of struct mtm_row in m.row_map is (re)sorted for either
	  * accessing the matrix by name or my offset.
	  */
	int                   offset;
	const char           *name;
	struct mtm_descriptor desc;
	MTM_ROW_PTR           data;
};

int  mtm_fetch_by_name( struct mtm_matrix *m, struct mtm_feature *f );
int  mtm_fetch_by_offset( struct mtm_matrix *m, struct mtm_feature *f );
const char *mtm_sclass_name( unsigned int );

extern const char *mtm_default_NA_regex;

int mtm_load_header( FILE *fp, struct mtm_matrix_header *header );
int mtm_load_matrix( FILE *fp, struct mtm_matrix *matrix, struct mtm_matrix_header *header );
#ifdef __cplusplus
}
#endif

#endif

