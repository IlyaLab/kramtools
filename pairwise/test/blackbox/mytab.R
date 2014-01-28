
## R's table function automagically does a lot of culling upon creation.
## I need to prevent that to simulate what pairwise does, thus my own
## table function is deliberately less "smart" than R's.

## WARNING: Resulting tables are only equivalent to pairwise' by virtue of
##          the fact that pairwise' pre-processor (prep.py) numbers 
##          categories alphabetically, the same as R.

.my.fill <- function(r,m) {
	if( ! any(is.na(r)) ) {
		rc <- as.integer(r);
		r  <- rc[[1]];
		c  <- rc[[2]];
		m[r,c] <- m[r,c] + 1;
		#tryCatch( m[r,c] <- m[r,c] + 1, error=function(e) {cat(sprintf("%d,%d\n",r,c))} );
	}
}

my.table <- function( df ) {

	t <- matrix( 0, nrow=length(levels(df[,1])), ncol=length(levels(df[,2])) );

	# Following removes factor-ish metadata so we have a simple matrix.
	# I was hitting trouble in which factors named "INF" were, in fact,
	# becoming Inf values (though I could not seem to reproduce this in
	# interactive mode).

	(apply( cbind( as.integer(df[,1]), as.integer(df[,2]) ), 1, .my.fill, t ));

	return(t);
}

