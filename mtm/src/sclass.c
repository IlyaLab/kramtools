
#include <ctype.h>
#include <assert.h>
#include "mtsclass.h"

static const char *STAT_CLASSES[ MTM_STATCLASS_COUNT ] = {
	"unknown",
	"boolean",
	"categorical",
	"ordinal",
	"continuous"
};


const char *mtm_sclass_name( unsigned int id ) {
	return id >= MTM_STATCLASS_COUNT ? "invalid" : STAT_CLASSES[ id ];
}

/**
  * These values CONSTRAIN the field type based on the specified
  * statistical class. They DETERMINE the field type when and only
  * when a single value is given (single bit is set).
  * This mapping embodies a set of conventions that are somewhat arbitrary,
  * but must be reflected in the implementation of _parseLine, particularly
  * the implementation of inference revision.
  *
  * Note that boolean and ordinal data both *could* be represented as
  * floating-point (e.g. 0.0/1.0), but these cases are unsupported now.
  * Keep it simple, stupid...
  *
  * We have encounted a real use-case in which boolean data was encoded
  * as {0,1} and the order mattered!
  */
static unsigned field_types[ MTM_STATCLASS_COUNT ] = {
	/* unknown     */ MTM_FIELD_TYPE_UNK,
	/* boolean     */ MTM_FIELD_TYPE_STR | MTM_FIELD_TYPE_INT,
	/* categorical */ MTM_FIELD_TYPE_STR,
	/* ordinal     */ MTM_FIELD_TYPE_INT,
	/* continuous  */ MTM_FIELD_TYPE_FLT
};

/**
  * Note that at runtime boolean, categorical and ordinal are encoded as
  * integers. Only continuous is encoded as floating-point.
  */

unsigned field_type_from_stat_class( unsigned id ) {
	return id >= MTM_STATCLASS_COUNT ? MTM_FIELD_TYPE_UNK : field_types[ id ];
}


/**
  * These are the conventions agreed on at ISB, by no means universal.
  * I'm making this conservative so that if the prefix does not really look
  * like the ISB convention described by the regex /[BCN]:/, then I assume 
  * that it's not really a type prefix and return UNKNOWN.
  */
int mtm_sclass_by_prefix( const char *token ) {

	if( ispunct( token[1] ) /* a colon at ISB, but allow any punct */ ) {
		switch( token[0] /* REQUIRE upper case */ ) {
		case 'B':
			return MTM_STATCLASS_BOOLEAN;
		case 'C':
		case 'F':
			return MTM_STATCLASS_CATEGORICAL;
		case 'D':
		case 'O':
			return MTM_STATCLASS_ORDINAL;
		case 'N':
			return MTM_STATCLASS_CONTINUOUS;
		}
	}
	return MTM_STATCLASS_UNKNOWN;
}

