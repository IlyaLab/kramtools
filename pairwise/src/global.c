
/**
  * These are not shared as widely as might be feared. They are
  * temporary expedients until the deferred reporting required for FDR 
  * and immediate reporting that is the default are fully rationalized.
  */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "stattest.h"
#include "global.h"
#include "format.h"

unsigned  arg_min_cell_count   = 5;
unsigned  arg_min_mixb_count   = 1;
unsigned  arg_min_sample_count = 2; // < 2 NEVER makes sense
double    arg_p_value          = 1.0;

FILE *g_fp_output = NULL;
FILE *g_fp_cache = NULL;

struct CovarsSummary g_summary;

row_t *g_rowmap = NULL;
size_t g_COLUMNS = 0;
const unsigned int *g_matrix_body = NULL;

void (*g_format)( 
		unsigned lf, 
		unsigned rf,
		unsigned status ) = format_tcga;

bool g_sigint_received = false;
