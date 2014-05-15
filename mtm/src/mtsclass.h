
#ifndef _mtsclass_h_
#define _mtsclass_h_

/**
  * I'm using bitfields below to accomodate the lack of mutual exclusion.
  * e.g. Categorical data may be boolean and/or ordered.
  */
#define MTM_STATCLASS_UNKNOWN      (0)
#define MTM_STATCLASS_BOOLEAN      (1)
#define MTM_STATCLASS_CATEGORICAL  (2)
#define MTM_STATCLASS_ORDINAL      (4)
#define MTM_STATCLASS_CONTINUOUS   (8)

const char *stat_class_name( unsigned bitfield );
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

