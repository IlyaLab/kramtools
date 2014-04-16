
#ifndef _stattest_h_
#define _stattest_h_

/**
  * Encapsulates the 3 elements common to all statistical tests and 
  * procedures: sample count, statistic, and probability of the statistic 
  * under NULL hypothesis assumptions.
  */
#define MAXLEN_STATRESULT_LOG (127)

struct Statistic {

	const char  *name;
	unsigned int sample_count;
	double       value;
	double       probability;

	/**
	  * This is for recording things like minimum cell count
	  * in a contingency table or number of ties in data.
	  * Content depends on .name.
	  */
	double       extra[ 4 ];

	/**
	  */
	char         log[ MAXLEN_STATRESULT_LOG+1 ];
};

#endif

