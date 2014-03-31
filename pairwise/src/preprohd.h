
#ifndef _preprohd_h_
#define _preprohd_h_

/**
  *  8 "PAIRWISE"
  *  4 version      == xxxx_MAJ_MIN_PAT
  *  4 header size  == 96
  *  4 rows
  *  4 DATA columns
  *  8 offset (from file start) to strings
  *  8 offset (from file start) to row map
  *  4 flags
  * 32 MD5 encoding )
  * 20 padding
  * == 96 bytes total
  */

struct preprocessed_input_header {
	char sig[8];                     // "PAIRWISE"

	unsigned int version;            // 0x00020000

	unsigned int header_size;        // 96

	unsigned int count_rows;         //
	unsigned int count_data_columns;

	unsigned long offset_names;
	unsigned long offset_map;

	unsigned int flags;

	char md5[32];

	char reserved[20];

} __attribute__((packed));

#define HEADER_FLAG_NAMES_SORTED (0x00000001)

#endif

