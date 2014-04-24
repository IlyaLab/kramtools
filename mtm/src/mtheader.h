
#ifndef _mtheader_h_
#define _mtheader_h_

/**
  * This defines the header for a saved mtmatrix.
  *  8 "PAIRWISE"
  *  4 endian       == 0x04030201
  *  4 version      == 0x_____MAJ_MIN_PAT
  *  4 header size  == 84
  *  4 datum_size;
  *  4 rows
  *  4 columns
  *  8 offset
  *  8 offset
  *  8 offset
  *  8 offset
  *  4 flags
  * 32 MD5 encoding )
  */

struct mtm_matrix_header {

	char sig[8];                  // "PAIRWISE"

	unsigned int endian;

	unsigned int version;         // 0x00020000

	unsigned int header_size;     // 84
	unsigned int datum_size;      // sizeof(float)

	unsigned int rows;
	unsigned int columns;

	unsigned long off_data;
	unsigned long off_descrip;
	unsigned long off_strings;
	unsigned long off_rowmap;

	unsigned int flags;           // none defined yet

	char md5[32];

} __attribute__((packed));

#define HEADER_FLAG_NAMES_SORTED (0x00000001)

#endif

