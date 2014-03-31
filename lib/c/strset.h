
#ifndef _hashtab_h_
#define _hashtab_h_

void * szs_create( unsigned int max, int dup );

#define SZS_ADDED          (+1)
#define SZS_PRESENT        ( 0)
#define SZS_ZERO_KEY       (-1)
#define SZS_TABLE_FULL     (-2)

int    szs_insert( void *, const char * str, int * );

#define SZS_RESIZE_FAILED  (-1)
#define SZS_OUT_OF_MEMORY  (-2)
int    szs_grow( void * );

#define SZS_NOT_FOUND      (-1)
int    szs_index( void *, const char * str );

int    szs_count( void * );
void   szs_clear( void * );
void   szs_destroy( void * );

#endif

