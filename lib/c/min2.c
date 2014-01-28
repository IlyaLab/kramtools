
#include "min2.h"

static void __attribute__((always_inline)) _swap( unsigned int* a, unsigned int* b ) {
	const unsigned int t = *a;
	*a = *b; 
	*b = t;
}


unsigned int min2ui( const unsigned int *arr, int n, unsigned int *second ) {

	unsigned int M = --n; // the (hypothetical) minimum
	unsigned int m = --n; // the (hypothetical) next-to-minimum
	if( arr[m] < arr[M] ) _swap( &m, &M );
	while( n-- > 0 ) {
		unsigned int K = arr[n];
		if( K < arr[m] ) {
			if( K < arr[M] ) {
				m = M;
				M = n;
			} else {
				m = n;
			}
		}
	}
	*second = m;
	return M;
}


#ifdef _UNITTEST_MIN2_
#include <stdio.h>
#include <stdlib.h>
int main( int argc, char *argv[] ) {

	int N 
		= argc > 1 
		? atoi( argv[1] ) 
		: 10 /* enough for interactive use */;
	int rd, i = 0;
	char *line = NULL;
	size_t n = 0;
	FILE *fp 
		= argc > 2
		? fopen( argv[2], "r" )
		: stdin;
	unsigned int *arr 
		= (unsigned int *)alloca( N*sizeof(unsigned int) );
	unsigned int min, penult_min;
	while( i < N && getline( &line, &n, fp ) > 0 ) {
		if( line[0] == '#' ) {
			fputs( line, stdout );
			continue;
		}
		rd = atoi( line );
		if( rd < 0 ) {
			fprintf( stderr, "Skipping %d. Unsigned only!\n", rd );
		} else
			arr[i++] = rd;
	}
	N = i; // ...how many we -actually- have, for remainder.
	free( line );
	fclose( fp );

	min = min2ui( arr, i, &penult_min );

	for(i = 0; i < N; i++ ) {
		const char *tag = "";
		if( i == min ) {
			tag = "minimum";
		} else
		if( i == penult_min ) {
			tag = "penultimate minimum";
		}
		printf( "%d %s\n", arr[i], tag );
	}

	return EXIT_SUCCESS;
}
#endif

