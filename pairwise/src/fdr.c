
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

#include "stattest.h"
#include "analysis.h"
#include "global.h"
#include "fdr.h"

struct FDRCacheRecord {
	double p;
	unsigned a,b;
} __attribute__((packed));
typedef struct FDRCacheRecord FDRCacheRecord_t;
typedef const FDRCacheRecord_t FDRCACHERECORD_T;

static int _cmp( const void *pvl, const void *pvr ) {

	FDRCACHERECORD_T *l = (FDRCACHERECORD_T*)pvl;
	FDRCACHERECORD_T *r = (FDRCACHERECORD_T*)pvr;

	if( l->p == r->p )
		return  0;
	else
		return (l->p < r->p) ? -1 : +1;
}


/**
  * Note that this method does NOT take a cache FILE* like fdr_postprocess
  * ONLY because it must follow a prototype that makes it interchangeable
  * with the _emit static in main. 
  * TODO: Revisit this.
  */
void fdr_aggregate( unsigned a, unsigned b, unsigned status ) {

	if( status == 0 && g_summary.common.P <= arg_p_value ) {
		struct FDRCacheRecord rec = {
			.p = g_summary.common.P,
			.a = a,
			.b = b
		};
		fwrite( &rec, sizeof(rec), 1, g_fp_cache );
	}
}


void fdr_postprocess( FILE *cache, double Q ) {

	const long size
		= ftell( cache );
	const unsigned COUNT
		= size / sizeof(struct FDRCacheRecord);
	struct FDRCacheRecord *prec, *sortbuf
		= calloc( size, sizeof(struct FDRCacheRecord) );
	const double RATIO
		= Q/COUNT;
	int n = 0;

	rewind( cache );
		fread( sortbuf, sizeof(struct FDRCacheRecord), COUNT, cache );
	rewind( cache );

	qsort( sortbuf, COUNT, sizeof(struct FDRCacheRecord), _cmp );

	prec = sortbuf;
	while( prec->p <= (n+1)*RATIO ) {

		const unsigned *pa = g_matrix_body + g_COLUMNS*prec->a;
		const unsigned *pb = g_matrix_body + g_COLUMNS*prec->b;

		clear_summary( &g_summary );

		const unsigned status = analysis_exec( 
				*pa, (const float*)(pa+1),
				*pb, (const float*)(pb+1),
				&g_summary );

		g_format( prec->a, prec->b, status );

		if( g_sigint_received ) {
			time_t now = time(NULL);
			fprintf( stderr, "# FDR postprocess interrupted @ %s", ctime(&now) );
			break;
		}
		prec += 1;
		n    += 1;
	}
}

