
-- This file contains working examples of Lua coroutines for specifying
-- feature pairs.
-- The executable passes the count of features (rows) to Lua coroutines
-- as the sole argument. Coroutines are free to use or ignore the value.

-- By default the executable looks for a function named pair_generator,
-- but you can tell it to look for a different function with the 
-- --coroutine/-c option.

-- Be aware that the executable expects 0-based offsets; the first row
-- in the matris is row 0. However, Lua uses 1-based indexing; the first
-- item in a Lua array is array[1].

-- This coroutine is an implementation of the exhaustive "all pairs"
-- case. 

function pair_generator( matrix_row_count )
	for i = 0,(matrix_row_count-1) do
		for j = i+1,(matrix_row_count-1) do
			coroutine.yield (i,j)
		end
	end
	-- Just "fall off the end"; no need to "return" anything else.
end

-- Lua has special syntax for iterating arrays similar to the items(),
-- keys() and values() methods on Python's dict type.

function cross_product( matrix_row_count )
	local L = {3,5,7,11}
	local R = {2,4,6,8}
	for i,l in ipairs(L) do
		for j,r in ipairs(R) do
			coroutine.yield( l, r )
		end
	end
end


-- You can define arbitrary globals that you then reference in 
-- functions/coroutines.

FIXED = 3

function one_versus_all( matrix_row_count )
	for i = 0,(matrix_row_count-1) do
		if i ~= FIXED then
			coroutine.yield (FIXED,i)
		end
	end
end


-- A 2-arg mathematical function (irrelevant to pairwise. part of
-- testing).

function f(x,y)
	return ( x^2 * math.sin(y) ) / (1-x)
end

