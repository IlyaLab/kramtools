#ifndef _featqual_h_
#define _featqual_h_

/**
  * This bitfield stores information about a feature's "quality,"
  * or lack thereof. 
  * It applies to both discrete and numeric, though the fields' 
  * interpretation is slightly different.
  */
struct feature_qualifier {

	unsigned int unused:9;

	/**
	  * For both categorical and numeric features this is the primary
	  * degeneracy flag. It indicates only one unique non-missing value was
	  * present in the feature--that is, the (non-missing scalars) of the 
	  * feature all had the same value.
	  *
	  * Alternatively, it can also indicated ALL scalars were missing which
	  * is itself a form of constant. The two possibilities are distinguished
	  * by looking at <missing>.
	  *
	  * In either case, if this bit is set the feature is probably unuseable.
	  */
	unsigned int constant:1;

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
	  * For BOTH feature types, a count of the missing values.
	  */
	unsigned int missing:16;
};

typedef struct feature_qualifier feature_qualifier_t;

#endif

