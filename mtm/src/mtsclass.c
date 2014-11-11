
#include <ctype.h>
#include <assert.h>
#include "mtsclass.h"

static const char *STAT_CLASSES[ 9 /* STAT_CLASS_ bitfield value */ ] = {
	"unk",
	"bin",
	"cat",
	"cat|bin",
	"ord",
	"ord|bin",
	"ord|cat",
	"ord|cat|bin",
	"con"
};


const char *stat_class_name( unsigned bitfield ) {
	return bitfield > 8
		? "invalid"
		: STAT_CLASSES[ bitfield ];
}

/**
  * These values CONSTRAIN the field type based on the specified
  * statistical class. They DETERMINE the field type when and only
  * when a single value is given (single bit is set).
  * This mapping embodies a set of conventions that are somewhat arbitrary,
  * but must be reflected in the implementation of _parseLine, particularly
  * the implementation of inference revision.
  */
static unsigned field_types[ 9 /* STAT_CLASS_ bitfield value */ ] = {
	MTM_FIELD_TYPE_UNK,                                // unknown
	MTM_FIELD_TYPE_STR,// | MTM_FIELD_TYPE_INT,           // boolean
	MTM_FIELD_TYPE_STR,// | MTM_FIELD_TYPE_INT,           // categorical
	MTM_FIELD_TYPE_STR,// | MTM_FIELD_TYPE_INT,           // boolean OR categorical
	MTM_FIELD_TYPE_INT, // implies floats are an error    ordinal
	MTM_FIELD_TYPE_INT, // implies floats are an error    ordinal OR boolean
	MTM_FIELD_TYPE_INT, // implies floats are an error    ordinal OR categorical
	MTM_FIELD_TYPE_INT, // implies floats are an error    ordinal OR boolean OR categorical
	MTM_FIELD_TYPE_FLT
};


unsigned field_type_from_stat_class( unsigned bitfield ) {
	return bitfield > 8
		? MTM_FIELD_TYPE_UNK
		: field_types[ bitfield ];
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

static const char *_DOMINANT_STATCLASS[5] = {
	"unknown",		// 0 => 0
	"boolean",		// 1 => 32-31 = 1
	"categorical",	// 2 => 32-30 = 2
	"ordinal",		// 4 => 32-29 = 3
	"continuous"	// 8 => 32-28 = 4
};

/**
  * This treats statistical classes as mutually exclusive
  */
const char *mtm_sclass_name( unsigned int bitflags ) {
	assert( (0xFFFFFFF0 & bitflags ) == 0 );
	return _DOMINANT_STATCLASS[ bitflags ? 32-__builtin_clz(bitflags) : 0 ]; 
}

