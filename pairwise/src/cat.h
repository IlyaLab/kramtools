
#ifndef _cat_h_
#define _cat_h_

/**
 * Encapsulates RAM management for counts and expectation matrices.
 * In particular, this is intended to be instantiated for the 
 * largest necessary table and reused (within the context of iteration).
 *
 * Provides several statistics and table validation methods.
 */
class CatCovars {

#if defined(_DEBUG) || defined(_UNITTEST_CAT_)
	friend std::ostream& operator<<( std::ostream& os, const CatCovars& ob );
#endif
	typedef unsigned int coord_t;
	typedef unsigned int count_t;
	typedef double       prob_t;

	/**
	 * Buffers are sized to accomodate these dimensions at construction
	 * and are not changeable thereafter.
	 */
	const unsigned int row_capacity, col_capacity;

	/**
	  * "decl(ared)_" because the bounds says nothing about the row/col 
	  * content--empty rows/columns may exist! I'm just driving this point
	  * home in these members' names because it has already burned me once!
	  */
	unsigned int decl_rows, decl_cols;

	/**
	  * The counts matrix is stored packed against the head of this
	  * (heap-allocated) buffer. In other words, though this buffer
	  * is -expected- to be larger than required (since the class 
	  * is designed to be iteratively re-used), only the first 
	  * (decl_rows*decl_cols) of the buffer are used for any given 
	  * matrix.
	  */
	count_t *counts;

	/**
	  * Following are derivative data...computed from counts.
	  */
	count_t *rmarg;
	count_t *cmarg;
	prob_t  *expect;

	count_t sample_count;
	prob_t  minimumExpected;

	/**
	 * This is used strictly by badCellExists
	 */
	count_t minimumAllowedCellCount;

	/**
	  * This class is stateful in the interest of efficiency
	  * because some methods share preliminary computations:
	  * 1) statistics require expectation, and 
	  * 2) expectation requires marginals
	  * (which may also be used incidentally for other purposes), so this 
	  * variable captures the progress towards complete calculations
	  * to enable implicit just-in-time calculation.
	  */
	enum State {
		Nothing = 0,
		Marginals,
		Expectation
	} calculated;

	bool badCellExists( unsigned int *off ) const;
	inline void coordsOfOffset( unsigned int off, unsigned int *r, unsigned int *c ) const;
	inline count_t& count( coord_t r, coord_t c );

#ifdef HAVE_MERGE_OPTION
	void mergeRows( unsigned int a, unsigned int b );
	void mergeCols( unsigned int a, unsigned int b );
#endif
	inline void cullRow( unsigned int r );
	inline void cullCol( unsigned int c );

	void recalcSampleCount();

public:
#ifdef HAVE_MERGE_OPTION
	enum {
		MergeSuccess,
		MergeFailure,
		MergeUnnecessary
	};
#endif
	CatCovars( unsigned int rcap, unsigned int ccap );
	~CatCovars();

	void setMinCellCount( unsigned n );

	/**
	 * Sets the bounds of the underlying matrix that will
	 * be used by subsequent code and zeroes JUST THAT PART
	 * of the whole table.
	 */
	void clear( coord_t r, coord_t c );
	bool complete();

	/**
	 * These accessors entirely (and losslessly(?)) encapsulate the pointer 
	 * arithmetic to reach into the matrix. It should not exist anywhere 
	 * else!
	 */
	inline void push( coord_t r, coord_t c );
	inline size_t size() const;
	inline count_t count( coord_t r, coord_t c ) const;
	inline prob_t expected( coord_t r, coord_t c ) const;
	inline bool is2x2() const;

	/**
	  * If after accumulation there are fewer than 2 rows or fewer than 2
	  * columns, the table is degenerate. Rows or columns with zero marginal
	  * sum are OK! For a contingency table that's bad but not -degenerate-.
	  */
	bool immediatelyDegenerate() const;
	bool eventuallyDegenerate() const;

	/**
	 * This carries out a fairly complex algorithm to eliminated
	 * cells with counts < minimumAllowedCellCount. 
	 * Generally speaking, either rows or columns are added together in a 
	 * way that minimizes the total alterations to the matrix and favors 
	 * keeping it square. See the code for the truly gory details.
	 */
#ifdef HAVE_MERGE_OPTION
	long merge( char *log, int loglen );
#endif
	/**
	 * Returns the indices of the minimal row and columns marginals.
	 * Can't be const because might require recalc of marginals.
	 */
	void minimalMarginals( unsigned int *rm, unsigned int *cm );

	/**
	 * Note that unlike merging cullBadCells cannot fail, so return is
	 * simply a count of total decimations (rows AND columns) that occurred.
	 */
	unsigned int cullBadCells( char *log, int buflen );

	/**
	  * Following two MUST be called in order before any of
	  * the statistics are computed. In other words,
	  *    !!!! THESE ARE NOT CALLED IMPLICITLY! !!!!
	  * (Statistical methods which -ought- to be const can't be unless
	  *  'expect', '[rc]marg' and friends are made mutable.)
	  */
	void calc_marginals();
	void calc_expectation();

	/**
	 * Available statistics
	 * Note: can't const g() and chi_square() since they can trigger JIT
	 * calculation of expectation...unless I make expect and co. mutable.
	 */
#ifdef HAVE_G_TEST
	int g( /* unfinished */ );
#endif
	int chi_square( struct CommonStats *, struct ChiSquareStats * );
	int fisher_exact( struct CommonStats *, struct FisherExactStats * ) const;
	double spearman_rho() const;

	enum {
		NoTransformation = 0,
		OnlyCull,
		OnlyMerge,
		CullThenMerge
	};
};

// Private inlines

void CatCovars::coordsOfOffset( unsigned int off, unsigned int *r, unsigned int *c ) const {
	*r = off / decl_cols;
	*c = off % decl_cols;
}

CatCovars::count_t& CatCovars::count( coord_t r, coord_t c ) {
	assert( r < decl_rows and c < decl_cols );
	return *(counts + r*decl_cols + c);
}

// Public inlines 

/**
 * Only this one is public in order that we may count the samples
 * as they're added...
 */
void CatCovars::push( count_t r, count_t c ) {

	assert( r < decl_rows );
	assert( c < decl_cols );

	count( r, c )++;

	// size and marginals are kept up to date...

	sample_count++;

	rmarg[r] += 1;
	cmarg[c] += 1;
}


size_t CatCovars::size() const {
	return sample_count;
}

CatCovars::count_t CatCovars::count( count_t r, count_t c ) const {
	assert( r < decl_rows and c < decl_cols );
	return *(counts + r*decl_cols + c);
}

CatCovars::prob_t CatCovars::expected( count_t r, count_t c ) const {
	assert( r < decl_rows and c < decl_cols );
	return *(expect + r*decl_cols + c);
}

bool CatCovars::is2x2() const {
	return (2 == decl_rows) and (2 == decl_cols);
}

#endif

