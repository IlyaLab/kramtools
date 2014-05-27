
#ifndef _cat_h_
#define _cat_h_

typedef unsigned int coord_t;

void cat_destroy( void *pv );
void *cat_create( unsigned int rcap, unsigned int ccap );

/**
 * Sets the bounds of the underlying matrix that will
 * be used by subsequent code and zeroes JUST THAT PART
 * of the whole table.
 */
void cat_clear( void *pv, coord_t r, coord_t c );

void cat_setMinCellCount( void *pv, unsigned n );

/**
 * These accessors entirely (and losslessly(?)) encapsulate the pointer 
 * arithmetic to reach into the matrix. It should not exist anywhere 
 * else!
 */
void   cat_push( void *pv, coord_t r, coord_t c );
size_t cat_size( void *pv );
bool   cat_complete( void *pv );
bool   cat_is2x2( void *pv );

/**
 * This function cannot fail, so return is
 * simply a count of total decimations (rows AND columns) that occurred.
 */
unsigned int cat_cullBadCells( void *pv, char *log, int buflen );

/**
 * Available statistics
 * Note: can't const g() and chi_square() since they can trigger JIT
 * calculation of expectation...unless I make expect and co. mutable.
 */
#ifdef HAVE_G_TEST
int cat_g( void *pv, struct Statistic * );
#endif
int cat_chi_square( void *pv, struct Statistic * );
int cat_fisher_exact( void *pv, struct Statistic * );
double cat_spearman_rho( void *pv );

#endif

