
#ifndef __rank_h__
#define __rank_h__

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Return values of rank_floats* functions contain these bit fields.
 */
#define RANK_STATUS_TIES  (0x00000001) // ...if ANY ties exist
#define RANK_STATUS_CONST (0x00000002) // ...if ALL values are ties.

void *rank_alloc( int n );
int   rank_floats( float *base, const unsigned int N, int normalize, void *buf );
int   rank_floats_strided( float *base, const unsigned int N, int stride, int normalize, void *buf );
void  rank_free( void * );

#ifdef __cplusplus
}
#endif

#endif
