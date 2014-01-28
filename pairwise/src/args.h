
#ifndef _args_h_
#define _args_h_

#ifdef __cplusplus
extern "C" {
#endif

/**
  * Command line options/arguments that are global.
  */
extern unsigned  arg_min_cell_count;
extern unsigned  arg_min_mixb_count;
extern unsigned  arg_min_sample_count; // < 2 NEVER makes sense
extern double    arg_p_value;

#ifdef __cplusplus
}
#endif

#endif
