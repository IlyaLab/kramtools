

ut_memmap : memmap.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_MEMMAP_  $^

ut_rank : rank.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_RANK_  $^

ut_ftest : fisher.c
	$(CC) -o $@ -O3 $^ -D_UNITTEST_FISHER_ -I$(HOME)/include -L$(HOME)/lib -lgsl -lgslcblas -lm

ut_min2 : min2.c
	$(CC) -o $@ -g -O0 -D_UNITTEST_MIN2_ $^ $(LDFLAGS)

ut_strset : strset.c contrib/fnv/hash_32.c 
	$(CC) -o $@ $(CFLAGS) -Icontrib -std=c99 -D_POSIX_C_SOURCE=200809L -D_UNIT_TEST_STRSET $^ -lm

ut_dsp : dsp.c
	$(CC) -o $@ -g -O0 $(CFLAGS) -D_UNITTEST_DSP_ $^

