
#ifndef _fp_h_
#define _fp_h_

#ifdef __cplusplus
extern "C" {
#endif

bool is_valid_fp( double val, const char *file, int line );

#ifdef __cplusplus
}
#endif

#ifdef HAVE_FP_PARANOIA
#define IS_VALID_FP(x) is_valid_fp((x),__FILE__,__LINE__)
#else
#define IS_VALID_FP(x) (true)
#endif

#endif

