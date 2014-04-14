
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

/**
  * "Sheila's format"
  * 1  -- feature A
  * 2  -- feature B
  * 3  -- analysis type (NN,CN,NC,CC,1X,X2)
  *      N: numeric
  *      C: categorical
  *      1X (in categorical pairs) feature B only had 1 (surviving) category
  *      X2 (in categorical pairs) feature A only had 1 (surviving) category
  * 4  -- Spearman correlation coefficient
  * 5  -- number of samples used in pairwise test
  * 6  -- -log10(p-value) (of whatever test was performed)
  * 7  -- count of unused samples in feature A
  * 8  -- -log10(p-value) of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature A
  * 9  -- count of unused samples in feature B
  * 10 -- -log10(p-value) of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature B
  */
static double __attribute__((always_inline)) _clamped( double prob ) {
	prob = -log10( prob );
	return isinf(prob) || (prob > 300 ) ? 300.0 : prob;
}


void format_tcga( 
		unsigned lf, 
		unsigned rf,
		unsigned status ) {

	double NLOGP[3];
	char rho[8];

	if( g_summary.spearman_rho == -2.0 )
		strcpy( rho, "NA" );
	else 
		sprintf( rho, "%+.3f", g_summary.spearman_rho );

	NLOGP[0] = _clamped( g_summary.common.P );
	NLOGP[1] = _clamped( g_summary.waste[0].common.P );
	NLOGP[2] = _clamped( g_summary.waste[1].common.P );

	if( g_rowmap ) {
		fprintf( g_fp_output, "%s\t%s\t",  g_rowmap[ lf ].name, g_rowmap[ rf ].name );
	} else {
		fprintf( g_fp_output, "%d\t%d\t",          lf      ,         rf       );
	}

	fprintf( g_fp_output, 
		"%s\t"
		"%s\t"
		"%d\t"
		"%.3f\t"
		"%d\t%.3f\t"
		"%d\t%.3f\t"
		"%s\n", 
		COVAR_TYPE_STR[ g_summary.kind ],
		rho, 
		g_summary.common.N, 
		NLOGP[0       ],
		g_summary.waste[0].unused, NLOGP[1],
		g_summary.waste[1].unused, NLOGP[2],
		g_summary.log );
}



/**
  * Standard format
  * 1  -- feature A
  * 2  -- feature B
  * 3  -- analysis_type(NN,CN,NC,CC,1X,X2):hypothesis_test
  *      N: numeric
  *      C: categorical
  *      1X (in categorical pairs) feature B only had 1 (surviving) category
  *      X2 (in categorical pairs) feature A only had 1 (surviving) category
  *      Hypothesis test is a short unambiguous string hyphenated onto the
  *      analysis type so that it's still one column (for now).
  * 4  -- Spearman correlation coefficient
  * 5  -- number of samples used in pairwise test
  * 6  -- p-value
  * 7  -- count of unused samples in feature A
  * 8  -- p-value of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature A
  * 9  -- count of unused samples in feature B
  * 10 -- p-value of Kruskal-Wallis test of difference 
  *       between used and unused parts of feature B
  * 11 -- log (of contingency table culling)
  */
void format_standard( 
		unsigned lf, 
		unsigned rf,
		unsigned status ) {

	if( g_rowmap ) {
		fprintf( g_fp_output, "%s\t%s\t",  g_rowmap[ lf ].name, g_rowmap[ rf ].name );
	} else {
		fprintf( g_fp_output, "%d\t%d\t",            lf       ,           rf        );
	}

	fprintf( g_fp_output, 
		"%s:%s\t"
		"%.3f\t"
		"%d\t"
		"%.3e\t"
		"%d\t%.3e\t"
		"%d\t%.3e\t"
		"%s\n", 
		COVAR_TYPE_STR[ g_summary.kind ], HypTestNames[ g_summary.test ],
		g_summary.spearman_rho, 
		g_summary.common.N, 
		g_summary.common.P,
		g_summary.waste[0].unused, g_summary.waste[0].common.P,
		g_summary.waste[1].unused, g_summary.waste[1].common.P,
		g_summary.log );
}


void format_abbreviated( 
		unsigned lf, 
		unsigned rf,
		unsigned status ) {

	const double rho = g_summary.spearman_rho;

	if( g_rowmap )
		fprintf( g_fp_output, "%s\t%s\t", g_rowmap[ lf ].name, g_rowmap[ rf ].name );
	else
		fprintf( g_fp_output, "%d\t%d\t",          lf       ,          rf        );

	fprintf( g_fp_output, 
		"%s\t"
		"%04x\t"
		"%d\t%.3e\t"
		"%d\t%.3e\t%d\t"
		"%d\t%.3e\t%d\t"
		"%.3f\t"
		"%s\n", 
		COVAR_TYPE_STR[ g_summary.kind ],
		status,
		         g_summary.common.N,          g_summary.common.P,
		g_summary.waste[0].common.N, g_summary.waste[0].common.P, g_summary.waste[0].unused,
		g_summary.waste[1].common.N, g_summary.waste[1].common.P, g_summary.waste[1].unused,
		isnan(rho) ? +nan("nan") : rho, 
		g_summary.log );
}

