
#ifndef _mix_h_
#define _mix_h_

#include <vector>
#include <utility>

/**
 * This class encapsulates statistics and helper methods for heterogeneous
 * data types meaning combinations of numeric/continuous and 
 * ordinal/categorical.
 */
class MixCovars {
#if defined(_DEBUG) || defined(_UNITTEST_MIX_)
	friend std::ostream& operator<<( std::ostream&, const MixCovars& );
#endif
	typedef std::pair<float,unsigned int> covars_t;
	typedef std::vector< covars_t > paired_samples_t;

	paired_samples_t samples;

	/**
	  * All category labels pushed into this accumulator must
	  * be in [0,CATEGORY_CAPACITY).
	 */ 
	const unsigned int CATEGORY_CAPACITY;
	const unsigned int MIN_BINARY_COUNT;

	/**
	  * This essentially bounds the allowed category LABELS. Any
	  * category label push'ed in must be in [0,expected_categories).
	  * ONLY THIS MEMBER SHOULD BE USED TO SIZE COUNTS ARRAYS!
	  */
	unsigned int expected_categories;

	/**
	  * This is an actual count of distinct labels OBSERVED by push.
	  */
	unsigned int observed_categories;

	/**
	  * A buffer allocated once in construction and reused for 
	  * all subsequent analyses. Reallocation is not currently
	  * supported.
	  */
	unsigned int *category_count;
	double       meanRank;     // necessarily of both covariates
	double       sum_sq_dev;   // used for Kruskal-Wallis and Spearman

	/**
	  * Following struct array is ONLY for Spearman rho calculation.
	  * Note that the number of non-empty categories after data entry
	  * says NOTHING about WHICH categories (indices) are non-zero.
	  * An earlier implementation was assuming category_count[0,1] > 0,
	  * BUT THIS NEED NOT BE TRUE! Thus, all this hack...
	  */
	struct {
		unsigned int index;
		unsigned int count;
		unsigned int meanRank;
	} edge[2];
	double sum_dev_prod; // used only for Spearman-rho

	/**
	  * Find the minimum and maximum non-empty categories.
	  */
	unsigned minCat() const;
	unsigned maxCat() const;
	unsigned int rank_sums( double *sums );

	MixCovars( const MixCovars& rhs );
	// copy construction illegal...to expensive, no need.

public:
	explicit MixCovars( unsigned int cap, unsigned int minb );
	~MixCovars();

	/**
	 * Semantics of this class follow those of vector, so
	 * clear removes elements without changing the (reserved) capacity.
	 */
	void reserve( unsigned int );
	void clear( unsigned cap );
	bool complete();

	inline void push( float num, unsigned int );
	inline size_t size() const;

	inline bool degenerate() const;
	inline bool categoricalIsBinary() const;

	/**
	 * Available statistics
	 * Reordering of data is a side affect of kruskal_wallis
	 * and mann_whitney.
	 */
	int kruskal_wallis( struct CommonStats *cs, struct KruskalWallisStats *ks );
#ifdef HAVE_MANN_WHITNEY
	int mann_whitney( struct CommonStats *cs, struct MannWhitneyStats *ms );
#endif
	/**
	  * This MUST be called AFTER kruskal_wallis or mann_whitney.
	  * ...at least as things are coded now.
	  */
	double spearman_rho() const;
};


void MixCovars::push( float num, unsigned int cat ) {

	assert( cat < expected_categories );

	if( 0 == category_count[ cat ] ) observed_categories++;

	category_count[ cat ] += 1;

	// Not updating edges in here because the number of conditionals
	// executed for sample counts > 32 exceeds the work to find the
	// edges post-sample accumulation.
	samples.push_back( covars_t( num, cat ) );
}

size_t MixCovars::size() const {
	return samples.size();
}

bool MixCovars::degenerate() const {
	return observed_categories < 2 or size() < 2;
}

bool MixCovars::categoricalIsBinary() const {
	return observed_categories == 2;
}

#endif

