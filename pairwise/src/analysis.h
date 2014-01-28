
#ifndef _analysis_h_
#define _analysis_h_

#ifdef __cplusplus
extern "C" {
#endif

int analysis_init( int columns );

#define FAIL_SAMPLES 0x00000001 // Too few samples after filtering.
#define FAIL_E_DEGEN 0x00000002 // Early degeneracy, before filtering.
#define FAIL_L_DEGEN 0x00000004 // Late degeneracy, after filtering.
#define FAIL_TOOMANY 0x00000008 // Too many categories
#define FAIL_MATH    0x00000010 // Something went south in the math.

/**
 * Returns non-zero on error and 0 otherwise.
 */
unsigned analysis_exec( unsigned int h1, const float *f1, unsigned int h2, const float *f2, 
		struct CovarsSummary * );

#ifdef __cplusplus
}
#endif

#endif

