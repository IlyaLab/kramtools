
#ifndef __fisher_h__
#define __fisher_h__

/**
  * Warning: NOT THREAD-SAFE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  */

#ifdef __cplusplus
extern "C" {
#endif

/**
  * Pre-reserves buffer space to preclude just-in-time allocation
  * (or RE-allocation).
  */
int fexact_reserve( unsigned int n );

/**
  * Forced immediate free'ing of private buffer.
  * Idempotent.
  */
void fexact_release( void );

/**
  * Calculates just the two-tail probability of a 2x2 counts table 
  * the same as or "more extreme" than that given using the hypergeometric
  * distribution.
  * This is intended to be used in the context of iteration since buffer
  * space is allocated just-in-time and reused by subsequent calls, and
  * freed only at the end of the process' execution.
  * Space can optionally be explicitly reserved and freed using the
  * companion APIs.
  */
double fexact_prob( unsigned int x, unsigned int m, unsigned int n, unsigned int k );

#ifdef __cplusplus
}
#endif
#endif

