
#ifndef _mtsclass_h_
#define _mtsclass_h_

/**
  * Statistical classes are mutually exclusive.
  * The statistical class of a feature guides both the parsing of its
  * data as well as the choice of statistical tests involving it.
  */
enum StatisticalClass {
	MTM_STATCLASS_UNKNOWN = 0,
	MTM_STATCLASS_BOOLEAN,     // Ordered IFF input data is integral.
	MTM_STATCLASS_CATEGORICAL, // ALWAYS UNORDERED.
	MTM_STATCLASS_ORDINAL,
	MTM_STATCLASS_CONTINUOUS,
	MTM_STATCLASS_COUNT
};

/**
  * Statistical class constrains the interpretation of textual data.
  */
unsigned field_type_from_stat_class( unsigned bitfield );

/**
  * Like statistical class, these are not mutually exclusive.
  * We use these identifiers to represent HOW WE WILL TREAT tokens,
  * not WHAT THE TOKENS ARE. Categorical data, for example, is typically
  * represented as string labels but may also be represented by small
  * integers.
  * Use of bit fields allows us to express that either representation is
  * acceptable.
  */
#define MTM_FIELD_TYPE_UNK  (0) // token contains...
#define MTM_FIELD_TYPE_STR  (1) // ...arbitrary alphabetic characters
#define MTM_FIELD_TYPE_INT  (2) // ...strictly digits of and the conventional
                                // prefix of some base with sign prefix
#define MTM_FIELD_TYPE_FLT  (4) // ...base-10 digits with decimal point and/or
                                // exponential notation with sign prefix
#endif

