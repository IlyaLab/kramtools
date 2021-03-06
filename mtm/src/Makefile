
############################################################################
# These variables are your opportunity to customize.

# Redefine this depending on where you want the results put.
INSTALL_LIB=$(HOME)/lib

# Define this if you want MD5 hashing of the input matrix
# built into the library.
HAVE_MD5=1

############################################################################

CFLAGS=-fPIC
ifdef DEBUG
CFLAGS+=-O0 -g -D_DEBUG -Wall
INFIX=-dbg
else
CFLAGS+=-O3 -DNDEBUG
LDFLAGS+=-Wl,-s
INFIX=
endif

LDFLAGS=-lm

BASENAME=mtm$(INFIX)
SHARED_LIB=lib$(BASENAME).so
STATIC_LIB=lib$(BASENAME).a
EXECUTABLES=ppm i2n

SRCLIB=../../lib/c
CONTRIB=$(SRCLIB)/contrib
CFLAGS+=-I$(SRCLIB)

CFLAGS+=-std=gnu99
CFLAGS+=-I$(SRCLIB) -I$(CONTRIB) 
ifdef HAVE_MD5
CFLAGS+=-DHAVE_MD5=1
endif

OBJECTS=parser.o \
	load.o \
	matrix.o \
	toktype.o \
	syspage.o \
	cardinality.o \
	sclass.o \
	specialc.o \
	feature.o \
	$(SRCLIB)/strset.o \
	$(CONTRIB)/fnv/hash_32.o

ifdef HAVE_MD5
OBJECTS+=$(CONTRIB)/md5/md5.o 
endif


all : $(EXECUTABLES) $(STATIC_LIB)

main.o    : mtmatrix.h mtheader.h syspage.h mterror.h
load.o    : mtmatrix.h mtheader.h syspage.h mterror.h
sclass.o  : mtsclass.h
toktype.o : mtsclass.h toktype.h
cardinality.o :
feature.o : feature.h mtmatrix.h toktype.h mterror.h mtsclass.h specialc.h \
	$(SRCLIB)/strset.h
parser.o  : syspage.h mtmatrix.h mtheader.h feature.h  mterror.h mtsclass.h specialc.h
syspage.o : syspage.h 
i2n.o     : syspage.h mtmatrix.h mtheader.h mterror.h
$(CONTRIB)/md5/md5.o : $(CONTRIB)/md5/md5.h

$(STATIC_LIB) : $(OBJECTS)
	ar rcs $@ $^ 

$(SHARED_LIB) : $(OBJECTS)
	$(CC) -shared -fPIC -o $@ $(CFLAGS) $^ $(LDFLAGS)

ppm : main.o $(STATIC_LIB)
	$(CC) -o $@ -static $(CFLAGS) $< -lm -L. -l$(BASENAME)

i2n : i2n.o $(STATIC_LIB)
	$(CC) -o $@ -static $(CFLAGS) $< -lm -L. -l$(BASENAME)

rmred : cull.o $(CONTRIB)/md5/md5.o $(STATIC_LIB)
	$(CC) -o $@ -static $(CFLAGS) cull.o $(CONTRIB)/md5/md5.o -lm -L. -l$(BASENAME) 

############################################################################
# Unit tests
############################################################################

ut-toktype : toktype.c
	$(CC) -o $@ -O0 -g -D_UNIT_TEST_TOKTYPE $^ -lm

ut-cardinality : cardinality.c
	$(CC) -o $@ $(CFLAGS) -D_UNIT_TEST_CARDINALITY $^ -lm

ut-feature : feature.c toktype.o cardinality.o sclass.o specialc.o \
		$(SRCLIB)/strset.o \
		$(CONTRIB)/fnv/hash_32.o
	$(CC) -o $@ $(CFLAGS) -std=gnu99 -DUNIT_TEST_FEATURE=1 -I$(SRCLIB) -I$(CONTRIB) $^ -lm

install : $(SHARED_LIB) $(STATIC_LIB)
	install $^ $(INSTALL_LIB)

#	install $^ $(INSTALL_LIB)

clean :
	rm -rf lib$(BASENAME)* $(EXECUTABLES) ut-* $$(find . ../../lib -name "*.o")

.PHONY : clean install

