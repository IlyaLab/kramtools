
#include <unistd.h> // for sysconf
#include "syspage.h"

/**
  * Returm the smallest size that is 
  * 1. larger than n AND
  * 2. a multiple of the system page size.
  */
size_t RT_PAGE_SIZE = 0;
size_t RT_PAGE_MASK = 0;

size_t page_aligned_ceiling( size_t n ) {

	if( RT_PAGE_SIZE == 0 ) {
		RT_PAGE_SIZE = sysconf(_SC_PAGESIZE);
		RT_PAGE_MASK = RT_PAGE_SIZE-1;
	}
	if( (n & RT_PAGE_MASK) != 0 )
		return (n & ~RT_PAGE_MASK) + RT_PAGE_SIZE;
	else
		return n; // ...a 1 in PAGE_SIZE chance!
}

