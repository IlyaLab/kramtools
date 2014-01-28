
#ifndef __mmap_h__
#define __mmap_h__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _mapped_file {
	char name[ FILENAME_MAX+1 ];
	int fd;
	void *mem;     // ptr to 1st byte of mapping
	size_t offset; // ...from beginning of file
	size_t len;    // ...in bytes of the mapping
} mapped_file_t;

/**
  * mmf_memmap is transactional. It succeeds or fails entirely, and
  * nothing is leaked or hanging if it fails!
  */
extern int mmf_memmap( mapped_file_t * /* in/out */, size_t filelen );
extern void mmf_munmap( mapped_file_t * );

#ifdef __cplusplus
}
#endif
#endif

