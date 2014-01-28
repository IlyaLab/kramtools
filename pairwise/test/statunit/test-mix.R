
A <- commandArgs( trailingOnly=TRUE );
K <- if(length(A) > 0 ) as.integer(A[[1]]) else 2;
N <- if(length(A) > 1 ) as.integer(A[[2]]) else 10;
F <- stdout();

# Following is designed to guarantee some ties.
# Note that for optimal-case testing these can be very efficiently
# filtered from the output by coreutils. In other words, it's easier 
# to postproc messy output to create a clean case than it is to make
# creation of messy vs. clean output a modal option in this code...
#
# Samples should contain ~%20 ties.

x <- data.frame(
	cat=as.integer( sample( 0:(K-1), N, replace=TRUE ) ),
	num=sample( runif(N/2)*100.0, N, replace=TRUE ) );

k   <- with( x, kruskal.test( num, cat ) );
rho <- cor( x$cat, x$num, method="spearman" );

if( TRUE )
	cat( sprintf( "# K=%f p-value=%f, spearman=%f\n", k$statistic, k$p.value, rho ), file=F );

# MannWhitney/Wilcox
#if( K == 2 ) {
#	u <- with( x, wilcox.test( num[cat==0], num[cat==1] ) );
#	cat( sprintf( "# U=%f p-value=%f, spearman=%f\n", u$statistic, u$p.value, rho ), file=F );
#}

write.table( x, F, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE );

