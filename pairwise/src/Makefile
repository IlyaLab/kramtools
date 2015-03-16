
############################################################################
# Product file names
############################################################################

# This is the sole place version is defined.
# The version.h header is then generated from these definitions.

EXECUTABLE_BASENAME=pairwise
VER_MAJ=2
VER_MIN=1
VER_FIX=2

ifdef RELCAN
RC=1
ifdef DEBUG
VER_TAG=rc$(RC)-dbg
else
VER_TAG=rc$(RC)-rel
endif
else
VER_TAG=
endif

#VERSION=-D_VER_MAJ=$(VER_MAJ)\
#	   	-D_VER_MIN=$(VER_MIN)\
#	   	-D_VER_FIX=$(VER_FIX)\
#		-D_VER_TAG=\"$(VER_TAG)\"

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

OBJECTS=main.o \
	fp.o \
	featpair.o \
	fixfmt.o \
	varfmt.o \
	analysis.c \
	cat.c \
	mix.c \
	num.c \
	usage_full.c \
	usage_short.c

SRCLIB=../../lib/c
CONTRIB=$(SRCLIB)/contrib
MD5DIR=$(CONTRIB)/md5

LIBOBJECTS=$(addprefix $(SRCLIB)/, dsp.o rank.o fisher.o min2.o)

EXECUTABLES=$(VERSIONED_EXECUTABLE)

############################################################################
# Compilation options
############################################################################

ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
MTM=mtm-dbg
else
CFLAGS+=-O3 -DNDEBUG
MTM=mtm
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
LINKTYPE=
else
LINKTYPE=-static
endif

############################################################################
# Primary

all : $(EXECUTABLES)

version.h : Makefile
	echo "#ifndef _version_h_"               > $@
	echo "#define _version_h_"              >> $@
	echo "#define VER_MAJOR ($(VER_MAJ))"   >> $@
	echo "#define VER_MINOR ($(VER_MIN))"   >> $@
	echo "#define VER_PATCH ($(VER_FIX))"   >> $@
	echo "#define VER_TAG (\"$(VER_TAG)\")" >> $@
	echo "#endif"                           >> $@

usage_full.c : usage_full.txt
	sed -f usage.sed $< > $@

usage_short.c : usage_short.txt
	sed -f usage.sed $< > $@

$(SRCLIB)/dsp.o :
$(SRCLIB)/rank.o : $(SRCLIB)/rank.h
$(SRCLIB)/fisher.o : $(SRCLIB)/fisher.h
$(SRCLIB)/min2.o : $(SRCLIB)/min2.h

analysis.o : stattest.h featpair.h analysis.h cat.h mix.h num.h args.h limits.h
bvr.o : bvr.h
cat.o : limits.h stattest.h cat.h fisher.h min2.h bvr.h
featpair.o : featpair.h
fixfmt.o : featpair.h stattest.h analysis.h fixfmt.h varfmt.h
fp.o : fp.h
main.o : featpair.h stattest.h analysis.h varfmt.h fixfmt.h limits.h version.h
mix.o : stattest.h mix.h bvr.h limits.h
num.o : rank.h stattest.h num.h
usage_full.o :
usage_short.o :
varfmt.o : stattest.h analysis.h varfmt.h featpair.h

pairwise : $(VERSIONED_EXECUTABLE)

$(VERSIONED_EXECUTABLE) : $(OBJECTS) $(LIBOBJECTS)
	$(CC) -o $@ $(LINKTYPE) $(CFLAGS) $^ $(LDFLAGS) -lgslcblas -lgsl -lm -l$(MTM) -ldl

# Following target will be eliminated away as soon as gratuitous C++ purged.

############################################################################
# Unit tests

UNITTESTS=ut_mix ut_cat ut_num ut_bvr

unittests : $(UNITTESTS)

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
	rm -rf $(EXECUTABLES) ut_* $(ARCHIVE_NAME).tar.gz usage_*.c *.o $$(find . ../../lib -name "*.o") version.h

.PHONY : clean archive

