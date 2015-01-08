
#ifndef _hashtab_h_
#define _hashtab_h_

typedef unsigned int (*string_hash_fx)(const char *str, unsigned int hval);

void * szs_create( unsigned int max, int dup, string_hash_fx fxn, unsigned int seed );

#define SZS_ADDED          (+1)
#define SZS_PRESENT        ( 0)
#define SZS_ZERO_KEY       (-1)
#define SZS_TABLE_FULL     (-2)

int    szs_insert( void *, const char * str, unsigned int * );


#define SZS_RESIZE_FAILED  (-1)
#define SZS_OUT_OF_MEMORY  (-2)
int    szs_grow( void * );


/**
  */
#define SZS_NOT_FOUND      (-1)
int    szs_tag( void *, const char * str );

int    szs_count( void * );
/**
  * Clears the content without freeing the array. Memory allocated for
  * individual strings -is- released.
  * This is expressly to reuse existing array for new content.
  */
void   szs_clear( void * );

/**
  * This frees -everything-.
  */
void   szs_destroy( void * );

#ifdef UNIT_AUTO_TEST
void   szs_dump( FILE * );
#endif

#endif

