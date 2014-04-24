
############################################################################
# These variables are your opportunity to customize.

# Redefine this depending on where you want the results put.
INSTALL_LIB=$(HOME)/lib

# Define this if you want MD5 hashing of the input matrix
# built into the library.
HAVE_MD5=

############################################################################

ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
else
CFLAGS+=-O3 -DNDEBUG
endif

SRCLIB=../../lib/c
CONTRIB=$(SRCLIB)/contrib
CFLAGS+=-I$(SRCLIB)

CFLAGS+=-std=gnu99
CFLAGS+=-I$(SRCLIB) -I$(CONTRIB) 
ifdef HAVE_MD5
CFLAGS+=-DHAVE_MD5
endif

SOURCES=mtparser.c\
	mtmatrix.c\
	toktype.c\
	syspage.c\
	cardinality.c\
	mtsclass.c\
	$(SRCLIB)/strset.c\
	$(CONTRIB)/fnv/hash_32.c

ifdef HAVE_MD5
SOURCES+=$(CONTRIB)/md5/md5.c 
endif

OBJECTS=$(addsuffix .o,$(basename $(SOURCES)))

all : mtproc install

libmtm.so : $(SOURCES)
	$(CC) -shared -fPIC -o $@ $(CFLAGS) $^ -lm

libmtm.a : $(OBJECTS)
	ar rcs $@ $^ 

mtproc : main.c syspage.c $(CONTRIB)/md5/md5.c libmtm.so
	$(CC) -o $@ $(CFLAGS) $^ -L$(INSTALL_LIB) -lm -lmtm

############################################################################
# Unit tests
############################################################################

ut-toktype : toktype.c
	$(CC) -o $@ -O0 -g -D_UNIT_TEST_TOKTYPE $^ -lm

ut-cardinality : cardinality.c
	$(CC) -o $@ $(CFLAGS) -D_UNIT_TEST_CARDINALITY $^ -lm

install : libmtm.so libmtm.a
	install $^ $(INSTALL_LIB)
	install $^ $(INSTALL_LIB)

clean :
	rm -rf libmtm.* mtproc ut-* $$(find . ../../lib -name "*.o")

.PHONY : clean install

