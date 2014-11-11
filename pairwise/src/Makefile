
############################################################################
# Product file names
############################################################################

EXECUTABLE_BASENAME=pairwise
VER_MAJ=2
VER_MIN=0
VER_FIX=0
ifdef DEBUG
VER_TAG=rc7-dbg
else
VER_TAG=rc7-rel
endif
VERSION=-D_VER_MAJ=$(VER_MAJ)\
	   	-D_VER_MIN=$(VER_MIN)\
	   	-D_VER_FIX=$(VER_FIX)\
		-D_VER_TAG=\"$(VER_TAG)\"

VERSIONED_EXECUTABLE=$(EXECUTABLE_BASENAME)-$(VER_MAJ).$(VER_MIN).$(VER_FIX)$(VER_TAG)

ARCHIVE_NAME=pairwise

############################################################################
# Dependencies
############################################################################

GSLINC=$(HOME)/include
LIBDIR=$(HOME)/lib

############################################################################
# Sources
############################################################################

SOURCES=main.c \
	fp.c \
	featpair.c\
	fixfmt.c \
	varfmt.c

SRCLIB=../../lib/c
CONTRIB=$(SRCLIB)/contrib
MD5DIR=$(CONTRIB)/md5

LIBSOURCES=dsp.c \
	rank.c \
	fisher.c \
	min2.c

SOURCES+=$(addprefix $(SRCLIB)/,$(LIBSOURCES))

$(warning $(SOURCES) )

EXECUTABLES=$(VERSIONED_EXECUTABLE)

UNITTESTS=ut_mix \
	ut_cat \
	ut_num \
	ut_bvr

############################################################################
# Compilation options
############################################################################

ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
else
CFLAGS+=-O3 -DNDEBUG
endif

CFLAGS+=-std=c99 
CFLAGS+=-I$(SRCLIB)
CFLAGS+=-I$(GSLINC)
CFLAGS+=-I../../mtm/src
CFLAGS+=-D_POSIX_C_SOURCE=200809L
# ...for getline declarations in main.c
CFLAGS+=-D_SVID_SOURCE
# ...for strdup declarations varfmt.c

LDFLAGS+=-L$(LIBDIR)

CFLAGS+=-DHAVE_LUA
LDFLAGS+=-llua

ifeq ($(shell uname),Darwin)
	STATIC=""
else
	STATIC="-static"
endif

############################################################################
# Primary

all : $(EXECUTABLES) $(UNITTESTS)

pairwise : $(VERSIONED_EXECUTABLE)

usage_full.c : usage_full.txt
	sed -f usage.sed $< > $@

usage_short.c : usage_short.txt
	sed -f usage.sed $< > $@

$(VERSIONED_EXECUTABLE) : $(SOURCES) analysis.c cat.c mix.c num.c usage_full.c usage_short.c
	$(CC) -o $@ $(STATIC) $(VERSION) $(CFLAGS) $^ $(LDFLAGS) -lgslcblas -lgsl -lm -lmtm

# Following target will be eliminated away as soon as gratuitous C++ purged.

unittests : $(UNITTESTS)

############################################################################
# Unit tests

ut_bvr : bvr.c
	$(CC) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_BVR_ $^

ut_mix : mix.c
	$(CC) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_MIX_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_cat : cat.c $(addprefix $(SRCLIB)/,fisher.c min2.c)
	$(CC) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_CAT_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_num : num.c fp.c $(addprefix $(SRCLIB)/,rank.c)
	$(CC) -o $@ -g -O0 $(CFLAGS) -D_UNITTEST_NUM_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_varfmt : varfmt.c
	$(CC) -o $@ -g -O0 $(CFLAGS) -D_UNIT_TEST_VARFMT $^ -lm

############################################################################
# General targets

archive : clean
	mkdir -p $(ARCHIVE_NAME)
	cp -r src doc $(ARCHIVE_NAME)
	tar cfz "$(ARCHIVE_NAME).tar.gz" $(ARCHIVE_NAME)
	rm -rf "./$(ARCHIVE_NAME)"

clean :
	rm -rf $(EXECUTABLES) ut_* $(ARCHIVE_NAME).tar.gz usage_*.c *.o $$(find . ../../lib -name "*.o")

.PHONY : clean archive

