
#ifndef _feature_h_
#define _feature_h_

/**
  * Holds buffers and state to support parsing feature data from a text file
  * and converting it to a runtime representation--that is, encoding it.
  */
struct feature {

	/**
	  * All features are expected to be the same length since they are
	  * rows in a table. This is the number of entries in buf below.
	  */
	int length;

	/**
	  * The same buffer is used for both encoding types (float/int) since
	  * they're the same size and never used simultaneously.
	  * This union merely obviates ugly casting.
	  */
	union {
		mtm_fp_t  *num;
		mtm_int_t *cat;
	} buf;

	int label_length;

	/**
	  * Every line's first field is a row label, not data.
	  */
	bool expect_row_labels;

	/**
	  * Any field which matches this regular expression is considered
	  * missing.
	  */
	const char *missing_data_regex;

	int (*interpret_prefix)(const char *);

	/***********************************************************************
	  * Following are strictly for categorical features.
	  */

	/**
	  * Maximum expected/supported cardinality of categorical features.
	  * This is primarily (only?) to appropriately size the hash table
	  * inside the string set.
	  */
	int max_cardinality;

	/**
	  * The string set used to map categorical features to integer values.
	  */
	void *category_labels;
};

int feature_count_fields( const char *pc, const char );

int feature_alloc_encode_state( struct feature * );
int feature_encode( char * const line,
		struct feature *f,
		struct mtm_descriptor *d );
void feature_free_encode_state( struct feature * );

#endif

