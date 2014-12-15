
-- This script generates the cross-product of (the indices of) two 
-- groups of rows:
-- 1) one group of rows at the start of the file
-- 2) and all rows following the start group.
--
-- For example:
-- N:gene1	ex1	ex2	ex3 ...
-- N:gene2	ex1	ex2	ex3 ...
-- ...
-- N:geneM	ex1	ex2	ex3 ...
-- C:mrkr1	gt1	gt2	gt3 ...
-- C:mrkr2	gt1	gt2	gt3 ...
-- ...
-- C:mrkrP	gt1	gt2	gt3 ...
--
-- NOTICE that you must edit the script below to specify the number of
-- rows constituting the first group because pairwise only passes the
-- total row count to the script.

function pair_generator( matrix_row_count )
	local first_group_size = 5
	for i=0,first_group_size-1 do
		for j=first_group_size,matrix_row_count-1 do
			coroutine.yield( i, j )
			--print( i, j )
		end
	end
end

