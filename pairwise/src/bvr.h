
/**
  * Binary Variables' Ranks
  */

#ifndef _bvr_h_
#define _bvr_h_

/**
  * When dealing with ranks of a binary variable, e.g.
  *            m         n    
  * value: 0 0 0 1 1 1 1 1
  *  rank: 1 2 3 4 5 6 7 8
  * ...if m is the number of zeros (3 above) and n is the total vector 
  * length (8 above) the mean of the ranks of the 1's is...
  */
inline float mean_rank_of_ties( unsigned m, unsigned n ) {
	return ( n*(n+1.0) - m*(m+1.0) )
		               /
		         ( 2.0 * (n-m) );
}

#endif

