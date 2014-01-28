
/**
  * This just lets me insert floating-point validity checks
  * in random places in the code.
  */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "fp.h"

bool is_valid_fp( double val, const char *file, int line ) {
	const char *ty = "";
	switch( fpclassify(val) ) {
	case FP_NAN:       ty = "NaN"; break;
	case FP_INFINITE:  ty = "Inf"; break;
	case FP_SUBNORMAL: ty = "Sub"; break;
	default:
		return true;
	}
	fprintf( stderr, 
		"warning: %s floating-point at %s:%d\n", ty, file, line );
	return false;
}

