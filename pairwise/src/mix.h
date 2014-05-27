
#ifndef _mix_h_
#define _mix_h_

void mix_destroy( void *pv );
void *mix_create( unsigned int sample_capacity, unsigned int category_capacity );

void mix_clear( void *pv, unsigned cap );
void mix_push(  void *pv, float num, unsigned int cat );
size_t mix_size( void *pv );
bool mix_complete( void *pv );


bool mix_degenerate( void *pv );
bool mix_categoricalIsBinary( void *pv );

/**
 * Available statistics
 * Reordering of data is a side affect of kruskal_wallis
 * and mann_whitney.
 */
int mix_kruskal_wallis( void *pv, struct Statistic * );
#ifdef HAVE_MANN_WHITNEY
int mix_mann_whitney( void *pv, struct Statistic * );
#endif

/**
  * This MUST be called AFTER kruskal_wallis or mann_whitney.
  * ...at least as things are coded now.
  */
double mix_spearman_rho( void *pv );

#endif

