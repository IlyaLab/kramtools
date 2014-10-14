
#include <alloca.h>

/**
  * Utility functions to count the number of distinct values in an array.
  */
static inline int _offset_of( 
		const unsigned int val, const unsigned int *buf, int len ) {
	while( len-- > 0 ) if( buf[len] == val ) break;
	return len;
}
#define _ARRAY_CONTAINS(v,b,l) (_offset_of((v),(b),(l)) >= 0)

/**
  * Return the cardinality of the integers in <buf> treating it as a set.
  * Values equal to the <NA> argument are entirely excluded from
  * consideration; they are entirely ignored.
  * If the cardinality is found to exceed <largest_of_interest>, then
  * <largest_of_interest>+1 is returned.
  * 0 is return iff the set (after excluding <NA>'s!) is empty.
  */
int cardinality( 
		const unsigned int *buf, int len, int largest_of_interest, const unsigned int NA ) {

	unsigned int *unique_values 
		= alloca( (largest_of_interest+1) * sizeof(unsigned int) );
	int cardinality = 0;

	for(int i = 0; i < len; i++ ) {
		if( buf[i] != NA ) {
			if( ! _ARRAY_CONTAINS( buf[i], unique_values, cardinality ) ) {
				unique_values[ cardinality++ ] = buf[i];
				if( cardinality > largest_of_interest ) break;
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
		printf( "cardinality( <array>, %d, %d, %d ) => %d\n",
			N, MAX, NA, cardinality( array, N, MAX, NA ) );
	} else
		errx( -1, "%s <NA value> <max> val1 [ val2 val3 ... ]", argv[0] );
	return EXIT_SUCCESS;
}
#endif

