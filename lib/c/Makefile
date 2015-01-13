
# iut: Interactive Unit Test
# aut: Automated Unit Test

iut_memmap : memmap.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_MEMMAP_  $^

iut_rank : rank.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_RANK_  $^

iut_ftest : fisher.c
	$(CC) -o $@ -O3 $^ -D_UNITTEST_FISHER_ -I$(HOME)/include -L$(HOME)/lib -lgsl -lgslcblas -lm

iut_min2 : min2.c
	$(CC) -o $@ -g -O0 -D_UNITTEST_MIN2_ $^ $(LDFLAGS)

iut_strset : strset.c contrib/fnv/hash_32.c 
	$(CC) -o $@ $(CFLAGS) -Icontrib -std=c99 -D_POSIX_C_SOURCE=200809L -D_UNIT_TEST_STRSET $^ -lm

aut_strset : strset.c contrib/fnv/hash_32.c 
	$(CC) -o $@ $(CFLAGS) -Icontrib -std=c99 -D_POSIX_C_SOURCE=200809L -DUNIT_AUTO_TEST $^ -lm

iut_dsp : dsp.c
	$(CC) -o $@ -g -O0 $(CFLAGS) -D_UNITTEST_DSP_ $^

