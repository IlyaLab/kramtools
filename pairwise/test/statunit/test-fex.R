
# Forcing output to stdout because that allows decimation of a particular
# pair (using grep and co.), guaranteeing a zero cell. This allows testing
# another scenario without adding code either here or to catcovar.cpp.

A <- commandArgs( trailingOnly=TRUE );
R <- if(length(A) > 0 ) as.integer(A[[1]]) else 3;
C <- if(length(A) > 1 ) as.integer(A[[2]]) else 3;
N <- if(length(A) > 2 ) as.integer(A[[3]]) else 16;
F <- stdout();

m <- matrix( 1+round(runif(R*C)*N), nrow=R, byrow=TRUE );

cat( sprintf( "# p-value=%f\n", fisher.test( m )$p.value ), file=F );
write.table( m, F, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE );

