
############################################################################
# Product file names
############################################################################

EXECUTABLE_BASENAME=pairwise
VER_MAJ=2
VER_MIN=0
VER_FIX=0
VERSION=-D_VER_MAJ=$(VER_MAJ) -D_VER_MIN=$(VER_MIN) -D_VER_FIX=$(VER_FIX)

VERSIONED_EXECUTABLE=$(EXECUTABLE_BASENAME)-$(VER_MAJ).$(VER_MIN).$(VER_FIX)

ARCHIVE_NAME=pairwise

############################################################################
# Dependencies
############################################################################

GSLINC=$(HOME)/include
GSLLIB=$(HOME)/lib

############################################################################
# Sources
############################################################################

SOURCES=main.c \
	fp.c \
	analysis.cpp \
	cat.cpp \
	mix.cpp \
	num.cpp

#	format.c

LIBDIR=../../lib/c
CONTRIB=$(LIBDIR)/contrib
MD5DIR=$(CONTRIB)/md5

LIBSOURCES=dsp.c \
	rank.c \
	fisher.c \
	min2.c

SOURCES+=$(addprefix $(LIBDIR)/,$(LIBSOURCES))

$(warning $(SOURCES) )

EXECUTABLES=$(VERSIONED_EXECUTABLE)

UNITTESTS=ut_mix \
	ut_cat \
	ut_num \
	ut_bvr \
	ut_dsp

############################################################################
# Compilation options
############################################################################

ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
else
CFLAGS+=-O3 -DNDEBUG
endif

CFLAGS+=-I$(LIBDIR)
CFLAGS+=-I$(GSLINC) -D_MAX_CATEGORIES=32 -DITER_SPEC_BUFFER_LEN=1024
LDFLAGS+=-L$(GSLLIB)

############################################################################
# Primary

all : $(EXECUTABLES) $(UNITTESTS)

pairwise : $(VERSIONED_EXECUTABLE)

usage_full.c : usage_full.txt
	sed -f usage.sed $< > $@

$(VERSIONED_EXECUTABLE) : $(SOURCES) usage_full.c
	$(CC) -o $@ -static $(VERSION) $(CFLAGS) $^ $(LDFLAGS) -lgslcblas -lgsl -lstdc++ -lm

unittests : $(UNITTESTS)

############################################################################
# Unit tests

ut_bvr : bvr.c
	$(CC) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_BVR_ $^

ut_mix : mix.cpp
	$(CXX) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_MIX_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_cat : cat.cpp $(addprefix $(LIBDIR)/,fisher.c min2.c)
	$(CXX) -o $@ -g -O0 -D_DEBUG -Wall $(CFLAGS) -D_UNITTEST_CAT_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_num : num.cpp fp.c $(addprefix $(LIBDIR)/,rank.c)
	$(CXX) -o $@ -g -O0 $(CFLAGS) -D_UNITTEST_NUM_ $^ $(LDFLAGS) -lgslcblas -lgsl -lm

ut_prep : prep.c $(LIBDIR)/strset.c $(CONTRIB)/md5/md5.c
	$(CC) -o $@ -std=gnu99 -O0 -g -I$(LIBDIR) -I$(CONTRIB) -DHAVE_CHECKSUM -D_UNIT_TEST_PREP $^ -lm

############################################################################
# General targets

archive : clean
	mkdir -p $(ARCHIVE_NAME)
	cp -r src doc $(ARCHIVE_NAME)
	tar cfz "$(ARCHIVE_NAME).tar.gz" $(ARCHIVE_NAME)
	rm -rf "./$(ARCHIVE_NAME)"

clean :
	rm -rf $(EXECUTABLES) ut_* $(ARCHIVE_NAME).tar.gz usage_full.c

.PHONY : clean archive

