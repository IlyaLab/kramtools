
#ifndef _iter_h_
#define _iter_h_

#ifdef __cplusplus
extern "C" {
#endif

int iter_parse( const char *outer, const char *inner );
bool iter_begin( int limit, int *offo, int *offi );
bool iter_next( int *offo, int *offi );

#ifdef _DEBUG
void iter_dump( const char *prefix, FILE *fp );
#endif

#ifdef __cplusplus
}
#endif
#endif
