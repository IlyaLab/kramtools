
#ifndef __min2_h__
#define __min2_h__

#ifdef __cplusplus
extern "C" {
#endif

/**
  * This quirky function finds the INDICES of the TWO LEAST
  * unsigned ints in an array (without the overhead of a sort).
  */
unsigned int min2ui( const unsigned int *arr, int n, unsigned int *second );

#ifdef __cplusplus
}
#endif
#endif

