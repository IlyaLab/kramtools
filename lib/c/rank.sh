
R --slave <<EOF
A <- commandArgs( trailingOnly=TRUE );
N <- if(length(A) > 0 ) as.integer(A[[1]]) else 10;
f <- if(length(A) > 1 )            A[[2]]  else ".1.tab";
x <- runif(10)*N
y <- sample( x, 2*N, replace=TRUE )
write.table( cbind( y, rank(y) ), f, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
EOF

# Then use the following bash to verify identical results.
# cut -f 1 .1.tab | ./ut_rank $(cat .1.tab | wc -l) > .2.txt
