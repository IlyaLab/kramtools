
############################################################################
# These variables are your opportunity to customize.

# Redefine this depending on where you want the results put.
INSTALL_LIB=$(HOME)/lib

# Define this if you want MD5 hashing of the input matrix
# built into the library.
HAVE_MD5=

############################################################################

BASENAME=libmtm
LDFLAGS=-lm

ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
INFIX=-dbg
else
CFLAGS+=-O3 -DNDEBUG
LDFLAGS+=-Wl,-s
INFIX=
endif

SHARED_LIB="$(BASENAME)$(INFIX).so"
STATIC_LIB="$(BASENAME)$(INFIX).a"

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

%.o : %.c
	$(CC) -c -o $@  -fPIC $(CFLAGS) $^

$(SHARED_LIB) : $(OBJECTS)
	$(CC) -shared -fPIC -o $@ $(CFLAGS) $^ $(LDFLAGS)

$(STATIC_LIB) : $(OBJECTS)
	ar rcs $@ $^ 

mtproc : main.c syspage.c $(CONTRIB)/md5/md5.c $(SHARED_LIB)
	$(CC) -o $@ $(CFLAGS) $^ -L$(INSTALL_LIB) -lm -lmtm

############################################################################
# Unit tests
############################################################################

ut-toktype : toktype.c
	$(CC) -o $@ -O0 -g -D_UNIT_TEST_TOKTYPE $^ -lm

ut-cardinality : cardinality.c
	$(CC) -o $@ $(CFLAGS) -D_UNIT_TEST_CARDINALITY $^ -lm

install : $(SHARED_LIB) $(STATIC_LIB)
	install $^ $(INSTALL_LIB)
	install $^ $(INSTALL_LIB)

clean :
	rm -rf $(BASENAME)* mtproc ut-* $$(find . ../../lib -name "*.o")

.PHONY : clean install

