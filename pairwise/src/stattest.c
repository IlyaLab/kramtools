
#include <string.h>
#include <math.h>
#include <stddef.h> // for offsetof
#include "stattest.h"

const char *COVAR_TYPE_STR[] = {
	"??",
	"NN",
	"CN",
	"NC",
	"CC",
	"?!",
	"!?"
};

const char *HypTestNames[HypothesisTestCount] = {
	"None",
	"FshrX",
	"ChiSq",
	"MannW",
	"KrskW",
	"Spear",
	"Pears"
};


void clear_summary( struct CovarsSummary *cs ) {

	// Initializing the p-values is a bit paranoid, but filtering requires
	// valid values in there. Initializing the log to a default non-empty
	// value so that _emit'ters don't have to check for NULL and column
	// counts are ALWAYS the same.

	memset( cs, 0, offsetof(struct CovarsSummary,log) );
	// ...don't bother zeroing entire log buffer!

	// All following is paranoia to insure initialization
	// within calculations.

	CLEAR_COMMON(          &(cs->common) );

	cs->waste[0].unused = -1;
	CLEAR_COMMON( &(cs->waste[0].common) );

	cs->waste[1].unused = -1;
	CLEAR_COMMON( &(cs->waste[1].common) );

	cs->spearman_rho = nan("nan");
	cs->log[0] = EOLOG;
	cs->log[1] = '\0' ;
}

