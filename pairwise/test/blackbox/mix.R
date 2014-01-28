
options(warn=-1); # to silence 
A <- commandArgs( trailingOnly=TRUE );
x <- read.table( A[[1]], head=T );
x <- x[ apply(x,1,function(r) {!any(is.na(r))}),];

# Mixed features will ALWAYS be named "N" and "C".
P <- tryCatch( with( x, kruskal.test(N,C)$p.value   ), error=function(e) {1.0} );

# Badly degenerate features may entail entire columns of NA.
# R seems to default assigning such features a class of logical which
# was causing trouble for droplevels. Better to just check explicitly...

rho <- if( ! all(is.na(x$C)) && ( length(levels(droplevels(x$C))) == 2 ) ) {
		with( x, cor( rank(C), N, method="spearman" ) )
	} else {
		NaN
	}

cat( sprintf( "%.6f\t%.6f\n", 
		if(is.nan(P) ) 1.0 else P, 
		if(is.na(rho)) NaN else rho ) )

