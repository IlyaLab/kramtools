
# Forcing output to stdout because that allows decimation of a particular
# pair (using grep and co.), guaranteeing a zero cell. This allows testing
# another scenario without adding code either here or to catcovar.cpp.

A <- commandArgs( trailingOnly=TRUE );
K <- if(length(A) > 0 ) as.integer(A[[1]]) else 3;
N <- if(length(A) > 1 ) as.integer(A[[2]]) else 32;
include.na <- (length(A) > 2 );

F <- stdout();

options( warn=-1);

# -1 means NA in this test.

x <- data.frame(
	cat.1=sample( if( include.na ) c(-1,0:(K-1)) else 0:(K-1), N, replace=TRUE ),
	cat.2=sample( if( include.na ) c(-1,0:(K-1)) else 0:(K-1), N, replace=TRUE ) );

S <- if( K==2 ) {
		with( x, fisher.test(cat.1, cat.2 ) ) 
	} else {
		with( x, chisq.test( cat.1, cat.2, correct=FALSE ) );
	}

Rho <- if( K==2 ) {
		with( x, cor( cat.1, cat.2, method="spearman") ) 
	} else {
		NA;
	}

cat( sprintf( "# p-value=%f, spearman=%f\n", S$p.value, Rho ), file=F );
write.table( x, F, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE );

