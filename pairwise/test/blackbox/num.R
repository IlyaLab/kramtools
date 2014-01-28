
# There are at last two ways this function can go off
# the rails:
# 1) NaN production
# 2) If one vector only contains 1 value, NA is returned.

options(warn=-1); # to silence "NaNs produced" warning.
A <- commandArgs( trailingOnly=TRUE );
x <- read.table( A[[1]],head=T )
x <- x[ apply(x,1,function(r) {!any(is.na(r))}),];

if( nrow(x) < 3 ) {
	cat( sprintf( "1.0\tnan\n" ) );
} else {
	res <- tryCatch( cor.test( x$N1, x$N2, exact=FALSE, method="spearman" ), 
			error=function(e) {list(p.value=1.0,estimate=NaN)} );

# Don't let NA every progogate outside...
	cat( sprintf( "%.3e\t%.3e\n", 
			if( is.na(res$p.value)  ) 1.0 else res$p.value, 
			if( is.na(res$estimate) ) NaN else res$estimate ) );
}

