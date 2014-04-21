
#ifndef _stattest_h_
#define _stattest_h_

/**
  * Encapsulates the 3 elements common to all statistical tests and 
  * procedures: sample count, statistic, and probability of the statistic 
  * under NULL hypothesis assumptions.
  */
#define MAXLEN_STATRESULT_LOG (127)
#define MAXNUM_EXTRA_VALUES   (4)

struct Statistic {

	const char  *name;
	unsigned int sample_count;
	double       value;
	double       probability;

	/**
	  * Following are for recording things like minimum cell count in a 
	  * contingency table or number of ties in data. Content depends on 
	  * Statistic.name.
	  */
	double       extra_value[ MAXNUM_EXTRA_VALUES ];
	const char * extra_name[  MAXNUM_EXTRA_VALUES ];

	/**
	  * The number of extra that are valid.
	  */
	int          extra_used;

	/**
	  */
	char         log[ MAXLEN_STATRESULT_LOG+1 ];
};

#endif

