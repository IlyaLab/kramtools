
#ifndef _mtheader_h_
#define _mtheader_h_

/**
  * This defines the header for a saved mtmatrix.
  * The file format is simple.
  * 1. There are five total sections including the header
  * 2. Each section is page aligned (with 0x00 padding between the end of
  *    one section and the start of the next).
  * 3. Sections are ordered according to the enumeration below.
  * 4. Only the rowmap section requires "relocation" upon loading.
  *    That is, it's string pointer are actually offsets all of which
  *    require the address of the ROWID section added to them.
  */

#define MTM_SIGNATURE ("MULTIMX")

enum SECTION {
	S_MATRIX = 0,	// a 2D array of mtm_int_t
	S_DESCRIPTOR,	// an array of struct mtm_descriptor
	S_ROWID,		// a packed sequence of NUL-terminated strings
	S_ROWMAP,		// an array of struct mtm_row_id, pointing into S_ROWID
	S_COUNT
};

struct section_descriptor {
	size_t size;   // actual size (not including tail padding)
	size_t offset; // from start of file
} __attribute__((packed));

struct mtm_matrix_header {

	char sig[8];                  // "MULTIMX\0"

	unsigned int endian;
	unsigned int version;         // 0x00020000
	unsigned int flags;
	unsigned int header_size;     // 84
	unsigned int datum_size;      // sizeof(float)
	unsigned int rows;
	unsigned int columns;

	struct section_descriptor section[ S_COUNT ];

	char md5[32];

} __attribute__((packed)); // 132 bytes

#define MTMHDR_ROW_LABELS_PRESENT (0x00000001)
#define MTMHDR_ROW_LABELS_LEXORD  (0x00000002)

#endif

