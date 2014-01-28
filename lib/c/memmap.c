
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h> // for O_RDONLY
#include <errno.h>
#include <string.h>
#include <unistd.h> // for close

#include "memmap.h"

/**
  * These routines just wrap the standard memory mapping calls in such
  * a way as to make them transactional. Either all steps of memory mapping
  * succeed or all steps are unwound for "completely" failure.
  */

int mmf_memmap( mapped_file_t *m, size_t filelen ) {

	/**
	 * Only stat the file if necessary to determine length.
	 */

	if( filelen == 0 ) {
		struct stat info;
		if( stat( m->name, &info ) ) {
			fprintf( stderr, 
				"error: failed stat'ing '%s'\n"
				"\t%s\n", 
				m->name, strerror(errno) );
			return errno;
		}
		if( info.st_size == 0 )
			return -1;
		filelen = info.st_size;
	}

	// Map the whole file unless a length was specified.

	if( m->len == 0 )
		m->len = filelen - m->offset;

	if( m->offset + m->len > filelen ) {
		fprintf( stderr, 
			"error: offset (%ld) + mapped len (%ld) exceeds file size (%ld)\n",
			m->offset, m->len, filelen );
		return -1;
	}

	m->fd = open( m->name, O_RDONLY );

	if( m->fd < 0 ) {

		fprintf( stderr, 
			"error: failed open'ing '%s'\n"
			"\t%s\n", 
			m->name, strerror(errno) );
		return errno;

	} else {

		m->mem = mmap( 0, m->len, PROT_READ, MAP_SHARED, m->fd, m->offset );
		if( MAP_FAILED == m->mem ) {
			fprintf( stderr, 
				"error: failed mmap'ing '%s'\n"
				"\t%s\n", 
				m->name, strerror(errno) );
			close( m->fd );
			m->fd = -1;
			return errno;
		}
	}

	return 0;
}


void mmf_munmap( mapped_file_t *m ) {

	if( m ) {
		if( m->mem ) {
			munmap( m->mem, m->len );
			m->mem = 0;
			m->len = 0;
		}
		if( ! ( m->fd < 0 ) ) {
			close( m->fd );
			m->fd = -1;
		}
	}
}

#ifdef _UNITTEST_MEMMAP_
#include <stdlib.h>
#include <stdio.h>
int main( int argc, char *argv[] ) {

	int err = 0;
	mapped_file_t mapping;
	unsigned char expect;
	size_t offset;

	memset( &mapping, 0, sizeof(mapping) );

	if( argc < 4 ) {
		fprintf( stderr, "%s <file> <byte offset> <expected byte>\n", argv[0] );
		exit(-1);
	}

	strcpy( mapping.name, argv[1] );
	offset = strtol( argv[2], NULL, strncmp( argv[2],"0x", 2 ) == 0 ? 16 : 10 );
	expect = strtol( argv[3], NULL, strncmp( argv[3],"0x", 2 ) == 0 ? 16 : 10 );

	err = mmf_memmap( &mapping, 0 );
	if( err == 0 ) {
		const unsigned char *pc = mapping.mem;
		fprintf( stdout, "byte @ 0x%0x == %02x, expected %02x.\n%s\n", 
				offset, pc[offset], expect,
			  	pc[offset] == expect ? "Success" : "Failed" );
		mmf_munmap( &mapping );
	} 

	return err ? EXIT_FAILURE : EXIT_SUCCESS;
}
#endif

