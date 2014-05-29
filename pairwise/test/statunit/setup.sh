R --slave -f test-cat.R --args 4 1000 > a
../../src/ut_cat 4 4 a 
R --slave -f test-num.R --args   1000 > b
../../src/ut_num 1000 b
R --slave -f test-mix.R --args 3 1000 > c
../../src/ut_mix 3 1000 c
