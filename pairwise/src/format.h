#ifndef _format_h_
#define _format_h_

void format_tcga( 
		FEATURE_PAIR_T *pair, 
		COVARIATEANALYSIS_T *covan,
		FILE *fp );

void format_abbreviated( 
		FEATURE_PAIR_T *pair, 
		COVARIATEANALYSIS_T *covan,
		FILE *fp );

void format_standard( 
		FEATURE_PAIR_T *pair, 
		COVARIATEANALYSIS_T *covan,
		FILE *fp );

#endif

