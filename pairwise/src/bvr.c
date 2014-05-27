
#ifdef _UNITTEST_BVR_
#include "bvr.h"
#include <stdio.h>
#include <stdlib.h>
int main( int argc, char *argv[] ) {
	const int M = argc > 1 ? atoi(argv[1]) : 1;
	const int N = argc > 2 ? atoi(argv[2]) : 1;
	fprintf( stdout, "mean_rank_of_ties( %d, %d ) => %.3f\n", M, N, MEAN_RANK_OF_TIES( M, N ) );
	return EXIT_SUCCESS;
}
#endif

