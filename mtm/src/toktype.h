
#ifndef _toktype_h_
#define _toktype_h_

extern const char *toktype_name[];

bool toktype_is_na_marker( const char *pc );
int  toktype_init( const char *na_expression );
/**
  * Always returns exactly one of the MTM_FIELD_TYPE_x bits.
  */
int  toktype_infer_narrowest_type( const char *pc, unsigned int *base );

#endif

