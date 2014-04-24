
-- This is an example of how an "all pairs" run could be achieved with
-- Lua script.
-- The pairwise executable automatically looks for a coroutine named 
-- "pair_generator" in any Lua script, so by so naming the coroutine
-- you don't need to specify a different function name.

-- You can define arbitrary globals that you then reference in 
-- functions/coroutines.

feature_count=4

function pair_generator( N )
	for i = 0,(N-1) do
		for j = i+1,(N-1) do
			coroutine.yield (i,j)
		end
	end
	-- Just "fall off the end"; no need to "return" anything else.
end

function one_versus_all()
	fixed = 3
	for i = 0,(feature_count-1) do
		if i ~= fixed then
			coroutine.yield (3,i)
		end
	end
end


-- A 2-arg mathematical function

function f(x,y)
	return ( x^2 * math.sin(y) ) / (1-x)
end

