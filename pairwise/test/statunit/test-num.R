
A <- commandArgs( trailingOnly=TRUE );
N <- if(length(A) > 0 ) as.integer(A[[1]]) else 10;
F <- stdout();

# Following is designed to guarantee some ties.
# Note that for optimal-case testing these can be very efficiently
# filtered from the output by coreutils. In other words, it's easier 
# to postproc messy output to create a clean case than it is to make
# creation of messy vs. clean output a modal option in this code...

S <- runif(3*N/4)*100.0;

x <- data.frame(
	l=sample( S, N, replace=TRUE ),
	r=sample( S, N, replace=TRUE ) );

cat( sprintf( "# rho(Pearson )=%f\n", with( x, cor( l, r, method="pearson" )  ) ), file=F );
cat( sprintf( "# rho(Spearman)=%f\n", with( x, cor( l, r, method="spearman" ) ) ), file=F );

write.table( x, F, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE );

