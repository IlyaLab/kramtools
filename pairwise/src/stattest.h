
#ifndef _stattest_h_
#define _stattest_h_

enum StatisticalTestName {
	UnknownHypTest = 0,
	FisherExact,
	ChiSquare,
	MannWhitney,
	KruskalWallis,
	Spearman,
	Pearson,
	HypothesisTestCount
};

/**
  * Encapsulates the 3 elements common to all statistical tests:
  * sample count, statistic, and probability of the statistic under NULL
  * assumptions.
  */
struct StatisticalTest {
	enum StatisticalTestName test;
	unsigned int N;
	double statistic;
	double P;
};

#endif

