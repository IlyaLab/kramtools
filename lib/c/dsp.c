
/**
  * Decimal Sequence Parser.
  * Any contiguous sequence of one or more decimal digits is recognized as
  * a decimal value, and ANYTHING other than decimal digits is treated as a
  * separator.
  */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

int get_base10_ints( FILE *fp, int *index, int n ) {

	int i = 0;
	int digits = 0;

	memset( index, 0, n*sizeof(int) );

	while( i < n ) {
		const int c
			= fgetc( fp );
		if( isdigit(c) ) {
			index[i] *= (10   );
			index[i] += (c-'0');
			++digits;
		} else
		if( digits > 0 ) {
			++i;
			digits = 0;
		} else
		if( EOF == c )
			return 0;
	}

	return i;
}

#ifdef _UNITTEST_DSP_

#include <stdlib.h>
#include <string.h>
#include <alloca.h>

int main( int argc, char *argv[] ) {
	const int N 
		= argc > 1 ? atoi(argv[1]) : 2;
	const size_t S
		= N * sizeof(int);
	int *i
		= alloca( S );
	memset( i, 0, S );
	while( get_base10_ints( stdin, i, N ) == N ) {
		int n = 0;
		printf( "%d", i[n++] );
		while( n < N ) printf( "\t%d", i[n++] );
		fputc( '\n', stdout );
		memset( i, 0, S );
	}
	return EXIT_SUCCESS;
}
#endif

