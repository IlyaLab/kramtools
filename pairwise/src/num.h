
#ifndef _num_h_
#define _num_h_

#include <vector>
#include <utility>

/**
 * This class encapsulates statistics and helper methods for continuous/
 * scalar data types.
 */
class NumCovars {

#if defined(_DEBUG) || defined(_UNITTEST_NUM_)
	friend std::ostream& operator<<( std::ostream&, const NumCovars& );
#endif
	std::vector<float> l, r;
	void *rank_scratch;

	NumCovars( const NumCovars& rhs );
	// copy construction illegal...to expensive, no need.

public:
	NumCovars();
	~NumCovars();

	/**
	 * Semantics of this class follow those of vector, so
	 * clear removes elements without changing the (reserved) capacity.
	 */
	void reserve( unsigned int );
	void clear();
	bool complete();
	inline void push( float , float );
	inline size_t size() const;

	/**
	 * Available statistics
	 */
#ifdef HAVE_SCALAR_PEARSON
	int pearson_correlation( struct CommonStats *, struct PearsonStats * );
#endif
	int spearman_correlation( struct CommonStats *, struct SpearmanStats * );
};


void NumCovars::push( float n1, float n2 ) {
	l.push_back( n1 );
	r.push_back( n2 );
}

size_t NumCovars::size() const {
	return l.size();
}

#endif

