
#ifndef _limits_h_
#define _limits_h_

/**
  * No categorical feature in the input may have more than this number of
  * distinct categories. Statistical tests on categorical data degrade 
  * rapidly as the category count increases unless sample count increases
  * faster; 32 is very reasonable for all sample counts we anticipate.
  *
  * This value is used throughout code to bound various memory allocations.
  */
#define MAX_ALLOWED_CATEGORY_COUNT (32)

#endif

