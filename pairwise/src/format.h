#ifndef _format_h_
#define _format_h_

void format_tcga( 
		const struct mtm_feature_pair *pair, 
		const struct CovariateAnalysis *covan,
		FILE *fp );

void format_abbreviated( 
		const struct mtm_feature_pair *pair, 
		const struct CovariateAnalysis *covan,
		FILE *fp );

void format_standard( 
		const struct mtm_feature_pair *pair, 
		const struct CovariateAnalysis *covan,
		FILE *fp );

#endif

