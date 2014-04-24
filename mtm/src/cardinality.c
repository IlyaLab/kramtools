
#include <alloca.h>

/**
  * Utility functions to count the number of distinct values in an array.
  */
static int __attribute__((always_inline)) _offset_of( 
		const unsigned int val, const unsigned int *buf, int len ) {
	while( len-- > 0 ) if( buf[len] == val ) break;
	return len;
}
#define _ARRAY_CONTAINS(v,b,l) (_offset_of((v),(b),(l)) >= 0)

int cardinality( 
		const unsigned int *buf, int len, int largest_of_interest, const unsigned int NA ) {

	unsigned int *values 
		= alloca( largest_of_interest * sizeof(unsigned int) );
	int cardinality = 0;

	for(int i = 0; i < len; i++ ) {
		if( buf[i] != NA ) {
			if( ! _ARRAY_CONTAINS( buf[i], values, cardinality ) ) {
				values[ cardinality++ ] = buf[i];
				if( cardinality == largest_of_interest ) break;
			}
		}
	}
	return cardinality; 
}

#ifdef _UNIT_TEST_CARDINALITY

#include <stdlib.h>
#include <stdio.h>
#include <err.h>

int main( int argc, char *argv[] ) {

	if( argc > 3 ) {
		const unsigned int NA 
			= strtol( argv[1], NULL, 0 );
		const unsigned int MAX
			= strtol( argv[2], NULL, 0 );
		const unsigned int N
			= argc - 3;
		unsigned int *array
			= alloca( N*sizeof(unsigned int) );
		for(int i = 0; i < N; i++ )
			array[i] = strtol( argv[3+i], NULL, 0 );
		printf( "max{ |array|, %d } => %d\n", MAX, cardinality( array, N, MAX, NA ) );
	} else
		errx( -1, "%s <NA value> <max> val1 [ val2 val3 ... ]", argv[0] );
	return EXIT_SUCCESS;
}
#endif

