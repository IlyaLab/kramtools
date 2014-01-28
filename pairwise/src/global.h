
#ifndef _global_h_
#define _global_h_

#include "args.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
  * Output streams.
  * g_fp_output is the primary output.
  *
  * The cache is used only when FDR control is activated.
  * It holds all potentially-significant results until post-processing
  * when significance is finally determined.
  */
extern FILE *g_fp_output;
extern FILE *g_fp_cache;

/**
 * ALL statistical results for a pair of covariates reside in g_summary.
 */
extern struct CovarsSummary g_summary;

/**
  * If a string table of row names was compiled into the binary matrix
  * it is accessed through this array.
  */
struct row {
	unsigned    offset;
	const char *name;
};
typedef struct row row_t;
typedef const row_t ROW_T;

extern row_t *g_rowmap;
extern size_t g_COLUMNS;
extern const unsigned int *g_matrix_body;

/**
  * This function determines the primary (textual) output format.
  */
extern void (*g_format)( 
		unsigned lf, 
		unsigned rf,
		unsigned status );

extern bool g_sigint_received;

#ifdef __cplusplus
}
#endif

#endif

