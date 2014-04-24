
#ifndef _mterror_h_
#define _mterror_h_

#define MTM_E_GENERIC    (-1)
/**
  * Unexpected type
  * Matrix structural
  * I/O failure
  */

#define MTM_E_IO          (-2)
#define MTM_E_SYS         (-2)

/**
  * A violation of manifest limits or assumptions built into the library
  * (e.g. NAN_AS_UINT) or implied by earlier arguments (e.g. category 
  * count).
  */

#define MTM_E_LIMITS         (-3)

/**
  * A field in the matrix was uninterpretable.
  */
#define MTM_E_FORMAT_FIELD   (-4)

/**
  * Unexpected column count or 
  */
#define MTM_E_FORMAT_MATRIX  (-5)

/**
  * The regex provided to define missing data failed to compile.
  */
#define MTM_E_INIT_REGEX  (-6)

/**
  * Caller passed in a NULL pointer.
  */
#define MTM_E_NULLPTR     (-7)

/**
  * Lookup of feature (by offset or name) failed.
  * If the lookup was by-offset, this implies out of bounds.
  */
#define MTM_E_NO_SUCH_FEATURE (-8)

#define MTM_E_NO_ROW_LABELS   (-9)

#endif

