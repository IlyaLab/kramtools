

ut_memmap : memmap.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_MEMMAP_  $^

ut_rank : rank.c
	$(CC) -o $@ $(CFLAGS) -D_UNITTEST_RANK_  $^

ut_ftest : fisher.c
	$(CC) -o $@ -O3 $^ -D_UNITTEST_FISHER_ -I$(HOME)/include -L$(HOME)/lib -lgsl -lgslcblas -lm

ut_min2 : min2.c
	$(CC) -o $@ -g -O0 -D_UNITTEST_MIN2_ $^ $(LDFLAGS)

