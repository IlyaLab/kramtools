
source('mytab.R')

options(warn=-1); # to silence "Chi-squared approximation may be incorrect" warning.
A <- commandArgs( trailingOnly=TRUE );
x <- read.table( A[[1]], head=T );
# NOT pre-culling NA's because we're simulate pairwise' approach.
#x <- x[ apply(x,1,function(r) {!any(is.na(r))}),];


exclude <- function( e, L ) {
	if( e==1 ) 
		2:L 
	else 
	if( e==L ) 
		1:(L-1)
	else
		c( 1:(e-1), (e+1):L)
}


cull.table <- function( t ) {

	ops <- character()
	m  <- unname(as.matrix(t));
	nr <- nrow(m); # these are updated progressively below
	nc <- ncol(m); # these are updated progressively below

	while( ((nr>2) || (nc>2)) && any(m<5) ) {

		# Determine the first row and first column ("first" of each 
		# merely to match the C code that's being tested.)
		w <- which( apply( m, 1, function(r) any(r<5) ) );

		if( length(w) > 0 ) { # any rows with bad cells?
			r <- min( w )
			c <- min( which( m[r,] < 5 ) ); # MUST exist given the above.

			# Cull whatever has the lowest marginal: row or column.

			if( sum(m[r,]) < sum(m[,c]) ) {
				# Cull the row, if possible, else a column
				if( nr > 2 ) {
					m <- m[ exclude(r,nr), ]
					nr <- nr - 1
					ops <- append( ops, sprintf("R%d",r-1) ); # -1 to convert to 0-based...
				} else {
					m <- m[, exclude(c,nc) ]
					nc <- nc - 1
					ops <- append( ops, sprintf("C%d",c-1) ); # ...used by pairwise.
				}
			} else {
				# Cull the col, if possible
				if( nc > 2 ) {
					m <- m[, exclude(c,nc) ]
					nc <- nc - 1
					ops <- append( ops, sprintf("C%d",c-1) );
				} else {
					m <- m[ exclude(r,nr), ]
					nr <- nr - 1
					ops <- append( ops, sprintf("R%d",r-1) );
				}
			}
		}
	}
	return( list(mat=m,log=ops) )
}

# Need to insure there -is- a table
# nrow(x) could be < 2, even 0.

DEGENERATE_OUTPUT <- "1.0\t-\n";
if( nrow(x) < 2 ) {
	cat( DEGENERATE_OUTPUT );
} else {

	t <- my.table(x);
	if( nrow(t) >=2 && ncol(t) >= 2 ) {
		cull.results <- cull.table( table(x) );

		l <- if( length( cull.results$log ) > 0 ) paste( cull.results$log, collapse="" ) else "-";
		P <- if( nrow(cull.results$mat) == 2 && ncol(cull.results$mat) == 2 ) {
			tryCatch( fisher.test( cull.results$mat )$p.value, error=function(e) {1.0} )
		} else {
			tryCatch( chisq.test( cull.results$mat )$p.value, error=function(e) {1.0} )
		}
		cat( sprintf( "%.6f\t%s\n", P, l ) );
	} else {
		cat( DEGENERATE_OUTPUT );
	}
}

