
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "mtmatrix.h"
#include "featpair.h"

int fetch_by_name( struct mtm_matrix *m, struct feature_pair *pair ) {
	int econd;
	econd= mtm_fetch_by_name( m, &(pair->l) );
	if( econd )
		return econd;
	econd= mtm_fetch_by_name( m, &(pair->r) );
	return econd;
}

int fetch_by_offset( struct mtm_matrix *m, struct feature_pair *pair ) {
	int econd;
	econd = mtm_fetch_by_offset( m, &(pair->l) );
	if( econd )
		return econd;
	econd = mtm_fetch_by_offset( m, &(pair->r) );
	return econd;
}

