
#ifndef _syspage_h_
#define _syspage_h_

#ifdef __cplusplus
extern "C" {
#endif

extern size_t RT_PAGE_SIZE;
extern size_t RT_PAGE_MASK;

size_t page_aligned_ceiling( size_t n );

#ifdef __cplusplus
}
#endif

#endif

