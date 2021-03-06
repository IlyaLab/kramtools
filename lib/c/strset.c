
/**
  * A SIMPLE string set implementation:
  */

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "strset.h"

struct entry {
	const char *str;
	unsigned tagval;
};

struct table {
	int capacity;
	int occupancy;
	unsigned int mask;
	int dup;

	string_hash_fx hash;
	unsigned int   seed;

	struct entry *array;
};

#define ENTRY_IS_OCCUPIED(t,p) ((t)->array[(p)].str != NULL)
#define ENTRY_IS_EMPTY(t,p)    ((t)->array[(p)].str == NULL)

static int _is_power_of_2( int v ) {
	return 0 == (v & (v-1));
}

static int _power_of_2_upper_bound( int v ) {
	int i = 0;
	while( (1<<i) < v && i < 31 ) i++;
	return (1<<i);
}


/**
  * Rehash all valid entries in the entry array. There is no reason why 
  * every entry that fit in the previous (half-sized) table should not go
  * in the new table.
  * Warning: Disable strdup'ing temporarily
  */
static int _rehash( struct entry *cur, int n, struct table *t ) {
	int error = 0;
	// Turn OFF string duplicating while rehashing...
	const int DUPING = t->dup;
	t->dup = 0;
	// ...because either:
	// 1. we weren't duplicating anyway, or
	// 2. we're rehashing strings already duplicated once that we own!
	// szs_insert doesn't know which of these is true; but we do here.
	for(int i = 0; i < n; i++ ) {
		if( cur[i].str != NULL ) {
			if( szs_insert( t, cur[i].str, NULL ) != SZS_ADDED ) {
				error = -1;
				break;
			}
		}
	}
	// WARNING: Don't return without resetting t->dup!
	t->dup = DUPING;
	return error;
}


/***************************************************************************
  * Public
  */

/**
  * Only RAM exhaustion can cause failure, so that is implied if return is 
  * NULL...no need for more error reporting.
  */

void * szs_create( unsigned int max, int dup, string_hash_fx fxn, unsigned int seed ) {

	struct table *table
		= calloc( 1, sizeof(struct table) );

	if( table ) {

		// occupancy already zeroed above.

		table->capacity
			= _power_of_2_upper_bound(max);

		assert( _is_power_of_2(table->capacity) );

		table->mask
			= (table->capacity-1);

		table->array
			= calloc( table->capacity, sizeof(struct entry) );

		table->dup
			= dup;

		table->hash
			= fxn;
		table->seed
			= seed;

		if( NULL == table->array ) {
			free( table );
			table = NULL;
		}
	}
	return table;
}


/**
  * Doubles the capacity of the table.
  * TODO: Currently indices are not stable--that is, after resizing strings
  * are mapped to different indices.
  */
int szs_grow( void *ht ) {

	struct table *t 
		= (struct table*)ht;

	struct entry *cur
		= t->array;

	t->array
		= calloc( t->capacity*2, sizeof(struct entry) );

	if( t->array ) {

		const int OCCUPANCY // Hold the current sizes in case of failure.
			= t->occupancy;
		const int CAPACITY
			= t->capacity;

		t->occupancy = 0;   // Update the table's vital stats...
		t->capacity *= 2;
		t->mask = (t->capacity-1);

		if( _rehash( cur, CAPACITY, t ) ) {
			// On failure reset everything to its original state.
			free( t->array );
			t->occupancy = OCCUPANCY;
			t->capacity  = CAPACITY;
			t->mask = (t->capacity-1);
			t->array = cur;
			return SZS_RESIZE_FAILED;
		}
		assert( OCCUPANCY == t->occupancy );
		free( cur );
		return 0;
	}
	return SZS_OUT_OF_MEMORY;
}


int szs_insert( void *ht, const char *str, unsigned int *tagval ) {

	struct table *t 
		= (struct table*)ht;

	assert( (t != NULL) && (str != NULL) );

	if( *str ) {

		const int K
			= t->hash( str, t->seed );
		const int IDEAL
			= K & t->mask;
		int pos
			= IDEAL;

		// Find either an empty slot or a slot with the same key.

		while( ENTRY_IS_OCCUPIED(t,pos) && strcmp( t->array[pos].str, str ) ) {
			pos = (pos + 1) % t->capacity;
			if( pos == IDEAL )
				break;
		}

		if( ENTRY_IS_EMPTY(t,pos) ) {
			t->array[pos].str   = t->dup ? strdup(str) : str;
			t->array[pos].tagval = t->occupancy++;
			if( tagval ) *tagval = t->array[pos].tagval;
			return SZS_ADDED;
		} else
		if( strcmp( t->array[pos].str, str ) == 0 ) {
			if( tagval ) *tagval = t->array[pos].tagval;
			return SZS_PRESENT;
		}
	} else
		return SZS_ZERO_KEY;

	return SZS_TABLE_FULL;
}


int szs_tag( void *ht, const char *str ) {

	struct table *t 
		= (struct table*)ht;

	assert( (t != NULL) && (str != NULL) );

	if( *str ) { // Zero should never be a key

		const int K
			= t->hash( str, t->seed );
		int remaining 
			= t->capacity;
		int pos
			= K & t->mask;
		while( strcmp( t->array[pos].str,str) && remaining-- > 0 ) 
			pos = (pos + 1) % t->capacity;
		if( strcmp( t->array[pos].str, str ) == 0 ) {
			return t->array[pos].tagval;
		}
	} else
		return SZS_ZERO_KEY;

	return SZS_NOT_FOUND;
}


int szs_count( void *ht ) {
	return ((struct table *)ht)->occupancy;
}


void   szs_clear( void *ht ) {
	struct table *t 
		= (struct table*)ht;
	if( t->dup ) {
		for(int i = 0; i < t->capacity; i++ ) {
			if( t->array[i].str )
				free( (void*)(t->array[i].str) );
		}
	}
	memset( t->array, 0, t->capacity*sizeof(struct entry) );
	t->occupancy = 0;
}


void szs_destroy( void *ht ) {

	struct table *t 
		= (struct table*)ht;

	assert( t != NULL );
	szs_clear( t );
	free( t->array );
	free( t );
}

#ifdef UNIT_AUTO_TEST

#include <stdio.h>
#include "fnv/fnv.h"

/**
  * Sort the struct entries of the hash table such that all entries with
  * valid strings precede all entries without strings (.str==NULL), and
  * secondarily by tag values.
  */
static int _cmp_entry( const void *pvl, const void *pvr ) {

	const struct entry *l = (const struct entry *)pvl;
	const struct entry *r = (const struct entry *)pvr;
	int result = 0;

	if( l->str != NULL ) {
		if( r->str != NULL ) {
			assert( l->tagval != r->tagval );
			result = (int)(l->tagval) - (int)(r->tagval);
		} else
			result = -1;
	} else
		result = r->str ? +1 : 0;

	return result;
}

int main( int argc, char *argv[] ) {

	int exit_status
		= 0;
	const int CAP
		= argc > 1 ? atoi( argv[1] ) : 10;
	const int DUP
		= argc > 2 ? atoi( argv[2] ) : 1;
	const int verbose
		= (int)( getenv("VERBOSE") != NULL );
	FILE *fp
		= argc > 3
		? fopen( argv[3], "r" )
		: stdin;

	char *line = NULL;
	size_t blen = 0;
	ssize_t llen;
	void *h;

	if( NULL == fp ) {
		fprintf( stderr, "no input\n" );
		exit(-1);
	}

	h = szs_create( CAP, DUP, fnv_32_str, FNV1_32_INIT );

	while( exit_status == 0 
			&& (llen = getline( &line, &blen, fp )) > 0 ) {
		int tagval;
		line[llen-1] = '\0'; // lop off the newline.
retry:
		switch( szs_insert( h, line, &tagval ) ) {
		case SZS_ADDED:
		case SZS_PRESENT:
			if( verbose )
				fprintf( stdout, "%s: OK(%d)\n", line, tagval );
			break;

		case SZS_ZERO_KEY:
			exit_status = 2;
			break;

		case SZS_TABLE_FULL:
			if( szs_grow( h ) )
				exit_status = 3;
			else
				goto retry;
			break;

		default:
			exit_status = 4;
		}
	}
	if( line )
		free( line );

	fclose( fp );

	if( exit_status == 0 ) {

		const int N
			= szs_count(h);
		struct table *t
			= (struct table*)h;
		int i;

		qsort( t->array, t->capacity, sizeof(struct entry), _cmp_entry );

		// Verify in a table with N entries,...
		// 1. they are tagged with values [0,n-1], and
		// 2. no other non-null entries exist in the table.

		for(i = 0; i < N; i++ ) {
			fprintf( stdout, "%s\n", t->array[i].str );
			if( i != t->array[i].tagval ) {
				exit_status = 5;
				goto EXIT;
			}
		}

		for(     ; i < t->capacity; i++ ) {
			if( t->array[i].str != NULL ) {
				exit_status = 6;
				goto EXIT;
			}
		}
	}

EXIT:
	szs_destroy( h );
	return exit_status;
}

#elif defined(_UNIT_TEST_STRSET)

#include <stdio.h>
#include <err.h>

#include "fnv/fnv.h"

int main( int argc, char *argv[] ) {

	if( argc > 1 ) {

		const int CAP
			= atoi( argv[1] );

		char *line = NULL;
		size_t blen = 0;
		ssize_t llen;
		void *h = szs_create( CAP, 1 /* duplicate */, fnv_32_str, FNV1_32_INIT );
		struct table *t
			= (struct table *)h;

		fprintf( stderr,
			"table capacity %d, mask %08x\n", t->capacity, t->mask );

		if( h ) {
			const char * k;
			unsigned int tagval, erv;
			while( (llen = getline( &line, &blen, stdin )) > 0 ) {
				if( line[0] == '?' ) {
					printf( "%d/%d\n", szs_count(h), t->capacity );
					continue;
				} else
				if( line[0] == '+' ) {
					if( szs_grow( h ) )
						errx( -1, "failed growing" );
					continue;
				} else
				if( line[0] == '-' ) {
					szs_clear( h );
					continue;
				}
				line[llen-1] = '\0'; // lop off the newline.
				switch( szs_insert( h, line, &tagval ) ) {
				case SZS_ADDED:
					printf( "+ %s => %d\n", line, tagval );
					break;
				case SZS_PRESENT:
					printf( "P %s => %d\n", line, tagval );
					break;
				case SZS_ZERO_KEY:
					printf( "! empty string\n" );
					break;
				case SZS_TABLE_FULL:
					printf( "! set full\n" );
					break;
				default:
					err( -1, "unexpected scanf count" );
				}
			}
			if( line )
				free( line );
			szs_destroy( h );
		} else
			err( -1, "failed creating hash table with capacity %d", CAP );
	} else
		errx( -1, "%s <capacity>", argv[0] );

	return EXIT_SUCCESS;
}
#endif

