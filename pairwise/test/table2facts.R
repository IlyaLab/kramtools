
# This script reads an integer matrix assumed to be a contingency
# table and translates the table into a 2-column table of (0-based)
# factor levels that would produce the same contingency table.
# i.e. It reverse transforms a contingency table into "raw" data.

A <- commandArgs( trailingOnly=TRUE );

i.file <- A[[1]];
o.file <- A[[2]];

x <- read.table( i.file, head=F)
sink( o.file )
for(r in 1:nrow(x) ) { 
	for(c in 1:ncol(x) ) { 
		for(n in 1:(x[r,c]) ) 
			cat( sprintf("%d\t%d\n", r-1, c-1 ) )
	}
}
sink()

