
#ifndef _num_h_
#define _num_h_

void   con_destroy( void *pv );
void  *con_create( unsigned int );
void   con_clear( void *pv );
void   con_push( void *pv, float n1, float n2 );
size_t con_size( void *pv );
bool   con_complete( void *pv );


#ifdef HAVE_SCALAR_PEARSON
int con_pearson_correlation( void *pv, struct Statistic * );
#endif

int con_spearman_correlation( void *pv, struct Statistic * );

#endif

