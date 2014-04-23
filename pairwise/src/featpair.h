
#ifndef _featpair_h_
#define _featpair_h_

/**
  * Since everything this application does involves pairs...
  * Some wrapper functions for the libmtm (single-feature) API.
  */

struct feature_pair {
	struct mtm_feature l, r; // left, right
};
typedef struct feature_pair feature_pair_t;
typedef const feature_pair_t FEATURE_PAIR_T;

extern int fetch_by_name( struct mtm_matrix *m, struct feature_pair *pair );
extern int fetch_by_offset( struct mtm_matrix *m, struct feature_pair *pair );

#endif

