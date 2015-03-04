
#ifndef _mtheader_h_
#define _mtheader_h_

/**
  * This defines the header for a binarized ("preprocessed") multi-type 
  * matrix saved to the filesystem. The file format is simple.
  * 1. There are four sections after the header
  * 2. Each section starts at an offset that is a multipleof the system's
  *    PAGE_SIZE (with 0x00 padding between the end of one section and the
  *    start of the next).
  * 3. Sections are ordered according to the enumeration below.
  * 4. Only the rowmap section requires "relocation" upon loading.
  *    That is, it's string pointers are actually offsets from the
  *    base of the ROWID section.
  *
  *    0x00000000 struct mtm_matrix_header
  *           ... 0x00 padding
  *    0x00001000 the data matrix[][]
  *           ... 0x00 padding
  *    0xUUUUU000 struct mtm_descriptor[]
  *           ... 0x00 padding
  *    0xVVVVV000 string1\0string2\0string3\0...
  *           ... 0x00 padding
  *    0xWWWWW000 struct mtm_row[]
  */

#define MTM_SIGNATURE ("MULTIMX")

/**
  * A preprocessed matrix stored in its binary form has S_COUNT sections.
  * The header is not considered a section because it's not loaded into the
  * single block of memory the way all sections are.
  */
enum SECTION {
	S_DATA = 0,	// a 2D array of mtm_int_t
	S_DESC,	    // an array of struct mtm_descriptor
	S_ROWID,	// a PACKED sequence of NUL-terminated strings
	S_ROWMAP,	// an array of struct mtm_row, pointing into S_ROWID
	S_COUNT
};

struct section_descriptor {
	size_t size;   // actual size (not including tail padding)
	size_t offset; // from start of file
} __attribute__((packed));

struct mtm_matrix_header {

	char sig[8];                  // "MULTIMX\0"

	unsigned int endian;
	unsigned int version;
	unsigned int flags;
	/**
	  * Size in bytes of this struct (in a file).
	  */
	unsigned int header_size;

	/**
	  * This is the span from the beginning of the S_DATA section (i.e. the
	  * first byte following the padding of the header *block*) to the last
	  * valid byte of the S_ROWMAP. In particular, the last section is not
	  * padded (on disk), so this size need NOT be a PAGE_SIZE multiple.
	  */
	size_t sizeof_rt_image;

	/**
	  * Size in bytes of one element of the matrix.
	  * Currently this is either sizeof(float) or sizeof(double).
	  */
	unsigned int sizeof_cell;
	unsigned int rows;
	unsigned int columns;

	struct section_descriptor section[ S_COUNT ];

	/**
	  * The digest of the original text input matrix (provenance).
	  */
	char md5[32+1];

} __attribute__((packed)); // 132 bytes

#define MTMHDR_ROW_LABELS_PRESENT (0x00000001)
#define MTMHDR_ROW_LABELS_LEXORD  (0x00000002)

#endif

