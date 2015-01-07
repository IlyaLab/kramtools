
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
  * actual data that is likely to be encountered/processed by this library.
  * In particular, I assume that it is "safe" to encode missing integral
  * data ("NAs") by the same -bit pattern- as (float)NaN (#defined below).
  *
  * For integer-encoded -categorical- data this is perfectly safe since
  * statistical considerations preclude categorical data from involving
  * more than (typically) 10's of categories. In other words, 0x7FC00000
  * is never likely to be needed to represent a category ("factor level").
  *
  * For ordinal data it's still safe, but maybe slightly less so.
  * Ordinal -should- rarely involve 2,143,289,344 levels (the decimal value
  * of 0x7FC00000), but I can't say that it can't.
  *
  * As a result of these considerations, I use the same -bit pattern- (NaN)
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
  * This bitfield stores information about a feature's "quality,"
  */
struct mtm_descriptor {

	unsigned int unused:7;

	/**
	  * This bit flags any row that, for any reason, should be not be used
	  * by the client. Currently, the only such reason is rows that are
	  * categorical with >32 categories.
	  */
	unsigned int ignore:1;

	/**
	  * For both integral and floating-point features this is the primary
	  * degeneracy flag. It indicates either:
	  * 1. All present fields' values were equal.
	  *    That is, only one unique value was present in the feature.
	  * 2. All fields' values were missing ("NA").
	  *    That is, the feauture is entirely empty of data.
	  *
	  * The two possibilities are distinguished by looking at <missing>.
	  * If missing == column count, it's case #2.
	  *
	  * In either case, if this bit is set the feature is probably unuseable.
	  */
	unsigned int constant:1;

	/**
	  * Set indicates integer data.
	  * This flag simply distinguishes between data encoded as float-point and
	  * data encoded as integer. Distinguishing between boolean, categorical
	  * and ordinal (all of which are integral) is another matter.
	  */
	unsigned int integral:1;

	/**
	  * This is both a flag and count. Non-zero indicates a categorical
	  * feature AND the number of distinct categories.
	  *
	  * Obviously a count of 1 implies <constant> above, however, both
	  * this and the preceding field will always be set independently
	  * and therefore always be consistent. In other words, you can
	  * rely on:
	  *            ( categories == 1 ) ==> ( constant == 1 )
	  */
	unsigned int categories:6;

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
	unsigned int missing:16;
};
typedef struct mtm_descriptor mtm_descriptor_t;
typedef const mtm_descriptor_t MTM_DESCRIPTOR_T;
typedef mtm_descriptor_t *mtm_descriptor_ptr;

#define MTM_MAX_MISSING_VALUES (65535)

////////////////////////////////////////////////////////////////////////////

struct mtm_row_id {
	/**
	  * This is the offset -in the matrix itself- (which is NEVER reordered).
	  */
	unsigned    rowoff;
	const char *string;
};
typedef struct mtm_row_id mtm_row_id_t;
typedef const mtm_row_id_t MTM_ROW_ID_T;

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

	size_t size;

	/**
	  * This is actually a mix of unsigned int's and float's,
	  * both of which are not coincidentally 4-bytes in size.
	  */
	mtm_int_t *data;

	struct mtm_descriptor *prop;

	/**
	  * This is only valid if the file was JIT preprocessed.
	  * Otherwise, the rownames reside with the main (memory-mapped)
	  * body of the file (and won't need to be indpendently freed).
	  */
	const char *row_names;

	/**
	  * This is always a free-standing array regardless of how the
	  * matrix was acquired.
	  */
	struct mtm_row_id *row_map;

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

	/**
	  * This is only filled in if the MD5 hash was compiled into libmtm.
	  * It is optional.
	  */
	char md5[33];
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
		struct mtm_matrix *m
		);

/**
  * The immediate motivation for this library is pairwise analysis of
  * features, so the following are convenience APIs for this use case.
  */
struct mtm_feature {
	/**
	  * This is the offset -in the matrix itself- (which is NEVER reordered).
	  * The array of struct mtm_row_id in m.row_map is (re)sorted for either
	  * accessing the matrix by name or my offset.
	  */
	int                   offset;
	const char           *name;
	struct mtm_descriptor prop;
	MTM_ROW_PTR           data;
};

int  mtm_fetch_by_name( struct mtm_matrix *m, struct mtm_feature *f );
int  mtm_fetch_by_offset( struct mtm_matrix *m, struct mtm_feature *f );
const char *mtm_sclass_name( unsigned int );

extern const char *mtm_default_NA_regex;

#ifdef __cplusplus
}
#endif

#endif

